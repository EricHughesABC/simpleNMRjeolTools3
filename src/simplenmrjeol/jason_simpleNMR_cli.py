"""
jason_simpleNMR_cli.py

The runnable program JASON's External Tools plugin invokes (via the
`simplenmr-jeol` console-script entry point installed by this package).
Resolves the input .jjh5 file, drives the dialog flow, builds the JSON
payload (via json_converter.jeolData), submits it to the simpleNMR
server, and opens the results viewer.

All the actual HDF5/RDKit data-extraction logic lives in
json_converter.py — this file is the orchestration/CLI layer only.
"""

from __future__ import annotations

from pathlib import Path
import sys
import os
import platform
# ── Fix a minimal PATH when launched by a GUI app (e.g. JASON) ──────────────
# Diagnostics showed JASON launches this script with a bare PATH
# (/usr/bin:/bin:/usr/sbin:/sbin) and no CONDA_PREFIX/CONDA_DEFAULT_ENV set,
# since GUI-launched processes on macOS don't inherit a login shell's
# environment the way a Terminal session does. sys.executable is unaffected
# (the OS resolves it directly), but PATH still matters for anything this
# process or its children shell out to. Fixed here, before any other
# imports (including json_converter's h5py/rdkit imports), by deriving the
# conda env's bin directory from sys.executable itself rather than
# hardcoding a machine-specific path.
_conda_bin = str(Path(sys.executable).parent)
_current_path = os.environ.get("PATH", "")
if _conda_bin not in _current_path.split(os.pathsep):
    os.environ["PATH"] = _conda_bin + os.pathsep + _current_path
os.environ.setdefault("CONDA_PREFIX", str(Path(sys.executable).parent.parent))

from typing import Optional
import json
from datetime import datetime
import fire

from qtpy.QtWidgets import QApplication, QMessageBox, QFileDialog

# simplenmr_builder is a hard dependency: fail loudly and immediately with
# a clear fix, rather than silently falling back to stale local duplicate
# code — see the 2026-08-28 project history for why that fallback pattern
# was removed everywhere. This one try/except covers both this file's own
# submission-layer imports AND json_converter's imports (importing
# .json_converter below triggers its simplenmr_builder imports too, so an
# ImportError from either surfaces here with the same clear message).
try:
    from .json_converter import jeolData
    from simplenmr_builder import ContractError
    from simplenmr_builder.gui.submission import (
        SubmissionOutcome,
        submit_to_server,
        check_user_registration,
        open_result_viewer_subprocess,
    )
except ImportError as e:
    print(
        "ERROR: simplenmr_builder[gui,viewer] is not installed in this "
        "environment. Install it with:\n"
        '    pip install -e "<path-to-simpleNMRbuilder>[gui,viewer]"\n'
        f"\nUnderlying import error: {e}"
    )
    sys.exit(1)

def init_qt_app():
    """Initialize Qt application if not already present"""
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    return app


def get_file_dialog():
    """Open a file dialog to select a .jjh5 NMR file"""

    file_path, _ = QFileDialog.getOpenFileName(
        None,
        "Select JEOL NMR File",
        str(Path.home()),  # Start in home directory
        "JEOL NMR Files (*.jjh5);;All Files (*)",
    )
    return Path(file_path) if file_path else None


def show_info_message(title: str, message: str, message_type=QMessageBox.Information):
    """Show an information message box."""
    msg_box = QMessageBox()
    msg_box.setIcon(message_type)
    msg_box.setWindowTitle(title)
    msg_box.setText(message)
    msg_box.exec_()


def validate_file(fn):
    """Validate that the file exists and has the correct extension"""
    if not fn:
        return False, "No file provided"

    fn = Path(fn)

    if not fn.exists():
        return False, f"File does not exist: {fn}"

    if fn.suffix.lower() != ".jjh5":
        return False, f"Invalid file extension. Expected .jjh5, got {fn.suffix}"

    return True, "File is valid"


def find_input_jjh5() -> Path | None:
    """
    Look for an auto-saved "input.jjh5" in the current working directory.

    When JASON's External Tools plugin launches this script, it sets the
    current working directory to a timestamped run folder (e.g.
    ".../External Tools/<toolName>/<timestamp>/") and saves the document
    there as "input.jjh5". Returns that path if it exists, else None.
    """
    candidate = Path.cwd() / "input.jjh5"
    if candidate.exists():
        return candidate
    return None


def commandline(fn=None):
    """
    Handle command line arguments for JEOL NMR file processing

    Resolution order:
        1. An explicit fn argument, if provided and valid.
        2. An auto-saved "input.jjh5" in the current working directory
           (the JASON External Tools plugin convention).
        3. A file picker dialog, as a last resort.

    Args:
        fn: Path to a JEOL .jjh5 NMR file (optional)
    """
    # Ensure we have a Qt application context
    # Convert to Path object if provided
    if fn:
        fn = Path(fn)
    else:
        # No explicit file given - check for JASON's auto-saved input.jjh5
        # in the current working directory before resorting to the dialog.
        auto_fn = find_input_jjh5()
        if auto_fn:
            print(f"Found input.jjh5 in current working directory: {auto_fn}")
            fn = auto_fn

    # Validate the provided file
    is_valid, message = validate_file(fn)

    if not is_valid:
        print(f"Invalid input ({message}): opening file dialog...")

        # Open file dialog
        fn = get_file_dialog()

        if not fn:
            # Show information dialog for no file selected
            show_info_message("No File Selected", "No file selected. Exiting.")
            sys.exit(1)

        # Validate the selected file
        is_valid, message = validate_file(fn)

        if not is_valid:
            # Show critical error dialog
            show_info_message(
                "File Validation Error",
                f"Selected file is invalid: {message}",
                QMessageBox.Critical,
            )
            sys.exit(1)

    return fn

def log_launch_diagnostics(launch_note: str = "") -> Optional[Path]:
    """Append process-launch environment details to a fixed diagnostics log.

    Written to a fixed location (~/simpleNMR_logs/) rather than anywhere
    derived from the .jjh5 path, since the whole point is comparing a
    JASON-launched run against a Terminal-launched run - the CWD/working
    directory differs between those two by definition, so a fixed path is
    the only way to reliably find both afterwards.

    Captures the full environment (sorted) plus a few fields called out
    specifically because they can confirm/rule out macOS App Sandbox
    inheritance from a GUI parent process (JASON) as the cause of the
    "mach_msg ... msg too large" renderer crash:
      - APP_SANDBOX_CONTAINER_ID: set only for sandboxed processes.
      - XPC_SERVICE_NAME / __CFBundleIdentifier: typically present when a
        process is launched by macOS's launchd/XPC machinery (i.e. from a
        GUI app like JASON) rather than an interactive shell.

    Never raises - a diagnostics failure must not block the main program.
    """
    try:
        log_dir = Path.home() / "simpleNMR_logs"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / "launch_diagnostics.log"

        timestamp = datetime.now().isoformat(timespec="seconds")
        highlight_keys = [
            "APP_SANDBOX_CONTAINER_ID",
            "XPC_SERVICE_NAME",
            "__CFBundleIdentifier",
            "PATH",
            "HOME",
            "TMPDIR",
            "DYLD_LIBRARY_PATH",
            "DYLD_FRAMEWORK_PATH",
            "DYLD_FALLBACK_LIBRARY_PATH",
            "CONDA_PREFIX",
            "CONDA_DEFAULT_ENV",
            "PYTHONPATH",
            "QT_API",
            "QT_PLUGIN_PATH",
        ]

        with open(log_path, "a", encoding="utf-8") as f:
            f.write(f"\n{'=' * 70}\n")
            f.write(f"Timestamp: {timestamp}\n")
            if launch_note:
                f.write(f"Note: {launch_note}\n")
            f.write(f"sys.executable: {sys.executable}\n")
            f.write(f"sys.argv: {sys.argv}\n")
            f.write(f"cwd: {os.getcwd()}\n")
            f.write(f"platform: {platform.platform()}\n")
            f.write(f"mac_ver: {platform.mac_ver()}\n")
            f.write(f"python: {sys.version}\n")

            f.write("\n-- Highlighted variables (sandbox/launch indicators) --\n")
            for key in highlight_keys:
                f.write(f"{key} = {os.environ.get(key, '<not set>')}\n")

            f.write("\n-- Full environment (sorted) --\n")
            for key in sorted(os.environ.keys()):
                f.write(f"{key} = {os.environ[key]}\n")

        print(f"Launch diagnostics written to: {log_path}")
        return log_path
    except Exception as e:
        print(f"Failed to write launch diagnostics log: {e}")
        return None


def main() -> int:
    """
    Entry point installed as the `simplenmr-jeol` console script (see
    pyproject.toml's [project.scripts]) — this is what JASON's External
    Tools plugin actually invokes. Also runnable directly via
    `python -m simplenmrjeol` (see __main__.py) or `python
    jason_simpleNMR_cli.py` for local testing.
    """

    log_launch_diagnostics()

    # local_remote = "http://127.0.0.1:5000"
    local_remote = "https://test-simplenmr.pythonanywhere.com"
    # local_remote = "https://simplenmr.pythonanywhere.com"

    # At the start of your main program
    # Return value intentionally unused — this call's only job is to
    # ensure a QApplication instance exists in this process before any
    # dialog (file picker, message box, the spectrum-assignment dialog
    # inside jeolData) tries to construct a widget. Results are shown via
    # open_result_viewer_subprocess (a separate process), not by calling
    # .exec() on this app directly.
    init_qt_app()

    print("SERVER ADDRESS:", local_remote)

    ml_address = f"{local_remote}/check_machine_learning"

    simpleNMR_address = f"{local_remote}/simpleMNOVA"

    jeol_fn = fire.Fire(commandline)

    # check if file exists
    if not jeol_fn.exists():
        print("File not found:", jeol_fn)
        show_info_message("No File Found", f"File not found: {jeol_fn}")
        sys.exit()

    jeol_data = jeolData(jeol_fn)
    # jeol_data.updateDatawithExptNames()

    jeol_data.choosePeakPickedSpectaforSimpleNMR()

    print("jeol_data\n", dir(jeol_data))

    if jeol_data.analysis_cancelled:
        print("Analysis was cancelled.")
        show_info_message(
            "Analysis Cancelled", "The analysis was cancelled by the user."
        )
        sys.exit()

    try:
        jeol_dict = jeol_data.createJsonDict()
    except ContractError as e:
        print(
            f"\nERROR: simplenmr_builder rejected this submission before it was "
            f"written: {e}\nThis means the server would also reject (or crash on) "
            f"this file — fix the underlying data issue rather than bypassing this "
            f"check.\n"
        )
        show_info_message(
            "Submission Rejected",
            f"The data could not be validated for submission:\n\n{e}",
        )
        sys.exit(1)

    print(dir(jeol_data))

    # in 1D nmrspectrum swap Delta1 and Delta2 peak values
    for key, spectrum in jeol_dict.items():
        if isinstance(spectrum, dict) and spectrum.get("datatype") == "nmrspectrum":
            if spectrum.get("type") == "1D":
                for pk_id, peak in spectrum.get("peaks", {}).get("data", {}).items():
                    peak["delta1"], peak["delta2"] = peak.get("delta2"), peak.get(
                        "delta1"
                    )

    jeol_json = json.dumps(jeol_dict, indent=4)

    # save JSON to file
    # create json file in same directory as jeol_fn using pathlib functions
    json_file_path = Path(
        jeol_fn.parent, jeol_fn.name.replace(jeol_fn.suffix, "_jeol_inputdata.json")
    )
    with open(json_file_path, "w") as json_file:
        json_file.write(jeol_json)

    for key in jeol_dict.keys():
        print(f"Key: {key}")

    # check if HSQC_0 exists in jeol_dict keys
    if "HSQC_0" not in jeol_dict.keys():
        show_info_message("Missing Data", "HSQC data set is required for analysis.")
        sys.exit()

    # Check user registration
    if not check_user_registration(ml_address):
        show_info_message("Registration Error", "Unable to verify registration.")
        # myGUIDATAwarn("Unable to verify registration. Please check your internet connection or contact support.")
    else:

        print("\n7. Submitting to simpleNMR Server...")
        submission = submit_to_server(jeol_dict, simpleNMR_address)

        if submission.outcome == SubmissionOutcome.SUCCESS:
            print("Analysis complete! Opening results viewer...")
            # Launched as a SEPARATE process (when simplenmr_builder is
            # available) rather than in-process — previously JASON's own
            # process ran MainWindow directly and called sys.exit(app.exec()),
            # which could require a manual Ctrl-C to actually terminate:
            # QtWebEngine can leave background threads alive even after
            # the Qt event loop returns. The shared subprocess entry
            # point (html_viewer.py's __main__) force-exits with
            # os._exit() the moment its window closes, which this
            # process doesn't get the benefit of when running MainWindow
            # in-process. wait=True blocks here until the viewer window
            # closes, so this script's own lifetime still visibly tracks
            # the viewer's, matching the previous UX, but now without the
            # hang. Confirmed 2026-08-27 fixing the equivalent issue for
            # the Bruker converter; applied here for consistency even
            # though it hadn't been separately reported for JEOL yet.
            open_result_viewer_subprocess(submission, wait=True)
            sys.exit(0)
        elif submission.outcome == SubmissionOutcome.DIAGNOSTIC_HTML:
            print(
                f"Server returned a diagnostic report (status {submission.status_code}). "
                "Opening it in the viewer..."
            )
            open_result_viewer_subprocess(submission, wait=True)
            sys.exit(0)
        else:
            # submit_to_server() has already shown the appropriate dialog
            # (network error / server error / registration required or
            # expired) - nothing further to do here except exit. The JSON
            # file was already saved locally regardless.
            print(f"Server submission did not succeed: {submission.outcome.value}")
            sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
