import fire
import sys
from pathlib import Path
from guidata import qapplication
from guidata.qthelpers import qt_app_context
from qtpy.QtWidgets import QFileDialog


def get_file_dialog():
    """Open a file dialog to select a .jjh5 NMR file"""
    with qt_app_context():
        file_path, _ = QFileDialog.getOpenFileName(
            None,
            "Select JEOL NMR File",
            str(Path.home()),  # Start in home directory
            "JEOL NMR Files (*.jjh5);;All Files (*)",
        )
        return Path(file_path) if file_path else None


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


def commandline(fn=None):
    """
    Handle command line arguments for JEOL NMR file processing

    Args:
        fn: Path to a JEOL .jjh5 NMR file (optional)
    """
    # Convert to Path object if provided
    if fn:
        fn = Path(fn)

    # Validate the provided file
    is_valid, message = validate_file(fn)

    if not is_valid:
        print(f"Invalid input: {message}")
        print("Opening file dialog...")

        # Open file dialog
        fn = get_file_dialog()

        if not fn:
            print("No file selected. Exiting.")
            sys.exit(1)

        # Validate the selected file
        is_valid, message = validate_file(fn)

        if not is_valid:
            print(f"Selected file is invalid: {message}")
            sys.exit(1)

    print(f"Processing JEOL NMR file: {fn}")

    # Add your NMR file processing logic here

    return fn


if __name__ == "__main__":
    print("Returned value:\n\t", fire.Fire(commandline))
