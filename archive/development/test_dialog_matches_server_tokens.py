"""
Regression test for the confirmed 2026-08-27 bug: chooseSpectra.py's
dialog offered dropdown labels ("Pureshift", "HSQCCLIPCOSY") that did not
match the server's literal NMREXPERIMENTS tokens ("H1_pureshift",
"HSQC_CLIPCOSY"), causing submitted spectra of those types to either
crash (KeyError in append_specdata's counter dict) or silently vanish
(never classify server-side, no error surfaced anywhere).

This test cross-checks the dialog's fixed dropdown list directly against
simplenmr_builder's NMREXPERIMENTS constant, so any future drift between
"what the dialog can produce" and "what the server actually recognizes"
fails loudly here instead of silently discarding real spectral data.

Requires simplenmr_builder to be importable (pip install -e it, or add
its path to PYTHONPATH) — see this repo's README for the dependency setup.
"""

import os
import sys
from pathlib import Path

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

try:
    from simplenmr_builder.constants import NMREXPERIMENTS
except ImportError:
    NMREXPERIMENTS = None

sys.path.insert(0, str(Path(__file__).parent.parent))  # allow importing chooseSpectra directly


@pytest.mark.skipif(NMREXPERIMENTS is None, reason="simplenmr_builder not installed")
def test_dialog_experiment_types_are_all_valid_nmrexperiments_tokens():
    from qtpy.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(sys.argv)
    from chooseSpectra import NMRExperimentDialog

    dialog_types = set(NMRExperimentDialog(experiment_names=[]).experiment_types) - {"SKIP"}
    server_tokens = set(NMREXPERIMENTS)

    unrecognized = dialog_types - server_tokens
    assert not unrecognized, (
        f"Dialog offers experiment type(s) the server's NMREXPERIMENTS list doesn't "
        f"recognize: {unrecognized}. A spectrum assigned to any of these would "
        f"never classify server-side and would be silently discarded. Server tokens: "
        f"{sorted(server_tokens)}"
    )
