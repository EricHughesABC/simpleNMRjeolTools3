# simpleNMRjeolTools_v6.py
#
# v6 changes:
#   - commandline() now looks for an auto-saved "input.jjh5" in the current
#     working directory before falling back to the file dialog. This supports
#     JASON's External Tools plugin, which sets the CWD to a timestamped run
#     folder (.../External Tools/<toolName>/<timestamp>/) containing
#     "input.jjh5" when it launches this script.

from pathlib import Path
import sys
import os
import platform
import uuid

# ── Fix a minimal PATH when launched by a GUI app (e.g. JASON) ──────────────
# Diagnostics showed JASON launches this script with a bare PATH
# (/usr/bin:/bin:/usr/sbin:/sbin) and no CONDA_PREFIX/CONDA_DEFAULT_ENV set,
# since GUI-launched processes on macOS don't inherit a login shell's
# environment the way a Terminal session does. sys.executable is unaffected
# (the OS resolves it directly), but PATH still matters for anything this
# process or its children shell out to. Fixed here, before any other
# imports, by deriving the conda env's bin directory from sys.executable
# itself rather than hardcoding a machine-specific path.
_conda_bin = str(Path(sys.executable).parent)
_current_path = os.environ.get("PATH", "")
if _conda_bin not in _current_path.split(os.pathsep):
    os.environ["PATH"] = _conda_bin + os.pathsep + _current_path
os.environ.setdefault("CONDA_PREFIX", str(Path(sys.executable).parent.parent))
from typing import Dict, List, Optional
from dataclasses import dataclass
from enum import Enum
from datetime import datetime
import requests
import webbrowser
import threading
import h5py
import numpy as np
import pandas as pd
import json
from rdkit import Chem
import fire

# from guidata.qthelpers import qt_app_context

from chooseSpectra import NMRExperimentDialog
from displayHTML_pyside import MainWindow

from qtpy.QtWidgets import (
    QProgressDialog,
    QApplication,
    QMessageBox,
    QFileDialog,
    QDialog,
)
from qtpy.QtCore import Qt


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


def extractC13predictionsJEOL(fn: Path, atoms_df: pd.DataFrame) -> pd.DataFrame:
    """
    Extracts predicted C13 NMR chemical shifts from a JEOL HDF5 file and combines them with atom information.

    Parameters
    ----------
    fn : Path
        Path to the JEOL HDF5 file containing NMR prediction data.
    atoms_df : pd.DataFrame
        DataFrame containing atom information, including atom indices and number of protons (nH).

    Returns
    -------
    pd.DataFrame
        DataFrame containing extracted C13 shift predictions with columns:
        - atom_idx: Index of the atom.
        - atomNumber: Atom number (1-based).
        - NumShifts: Number of predicted shifts for the atom.
        - Shifts: Array of predicted C13 shifts.
        - CalcMethod: Calculation method used for prediction.
        - numProtons: Number of protons (from atoms_df) associated with each atom.

    Notes
    -----
    - If the required HDF5 keys are missing, a KeyError is caught and printed.
    - The function assumes a specific structure for the JEOL HDF5 file.
    """

    C13SHIFTS = 0

    c13shifts_df = pd.DataFrame(
        columns=["atom_idx", "atomNumber", "NumShifts", "Shifts", "CalcMethod"]
    )

    with h5py.File(fn, "r") as fp:
        try:
            c13shifts = fp[
                f"JasonDocument/Molecules/Molecules/0/NMRData/Spectra/SpectraList/{C13SHIFTS}/Shifts"
            ]

            print("C13 Shifts\n\t", list(c13shifts))
            c13shifts_data = []
            for shift in c13shifts:
                print(shift)
                numC13shifts = len(c13shifts[shift].attrs.get("Value", None))
                calcMethods = c13shifts[shift].attrs.get("Value.Method", None)
                # add row to dataframe
                c13_row = {
                    "atom_idx": int(c13shifts[shift].attrs.get("Nums", None)[0]),
                    "atomNumber": int(c13shifts[shift].attrs.get("Nums", None)[0]) + 1,
                    "NumShifts": numC13shifts,
                    "Shifts": c13shifts[shift].attrs.get("Value", None),
                    "CalcMethod": calcMethods,
                }
                c13shifts_data.append(c13_row)

            c13shifts_df = pd.DataFrame(c13shifts_data)
            c13shifts_df["numProtons"] = 0  # Initialize numProtons column

            for row in atoms_df.itertuples():
                atom_idx = row.atom_idx
                numProtons = int(row.nH)
                # replace numProtons in c13shifts_df using atom_idx

                idx = c13shifts_df.index[c13shifts_df["atom_idx"] == atom_idx].tolist()
                if idx:
                    c13shifts_df.loc[idx, "numProtons"] = numProtons
            #  convert numProtons to integer
            c13shifts_df["numProtons"] = c13shifts_df["numProtons"].astype(int)
        except KeyError as e:
            print(f"KeyError: {e}")

    return c13shifts_df


def create_c13predictions(c13shifts_df: pd.DataFrame) -> dict:
    """
    Generate a dictionary containing C-13 NMR shift predictions from a DataFrame.
    Args:
        c13shifts_df (pd.DataFrame): A pandas DataFrame with columns 'atom_idx', 'atomNumber',
            'numProtons', and 'Shifts'. The 'Shifts' column should contain a list of shift values.
    Returns:
        dict: A dictionary with the following structure:
            {
                "count": <number of atoms>,
                "data": {
                    <atom_idx>: {
                        "atom_idx": <int>,
                        "atomNumber": <int>,
                        "numProtons": <int>,
                        "ppm": <float>  # last value from 'Shifts' list
                    },
                    ...
    """

    c13predictions = {
        "datatype": "c13predictions",
        "count": len(c13shifts_df),
        "data": {},
    }

    for idx, row in c13shifts_df.iterrows():
        c13predictions["data"][row["atom_idx"]] = {
            "atom_idx": row["atom_idx"],
            "atomNumber": row["atomNumber"],
            "numProtons": row["numProtons"],
            "ppm": float(row["Shifts"][-1]),
        }

    return c13predictions


def get_host_name() -> dict:
    """
    Retrieves the MAC-based host identifier and returns it in a structured dictionary.

    Returns:
        dict: A dictionary containing the MAC-based host identifier under the key "data",
              along with metadata fields "datatype" and "count".
              Example:
              {
                  "datatype": "hostname",
                  "count": 1,
                  "data": {
                      "0": "<mac_based_id>"
                  }
              }
    """

    mac_based_id = hex(uuid.getnode())

    hostname = {"datatype": "hostname", "count": 1, "data": {"0": mac_based_id}}

    return hostname


# function to find datasets with peaks
def find_datasets_with_peaks(file_path: Path) -> list[str]:
    """
    Scans an HDF5 file for NMR experiment datasets containing peaks.

    Args:
        file_path (Path): Path to the HDF5 file to be scanned.

    Returns:
        dict: A dictionary where keys are experiment IDs with detected peaks,
              and values are dictionaries containing placeholders for
              'spec_info', 'peaks', 'integrals', and 'multiplets'.

    Side Effects:
        Prints the number of peaks found for each experiment, or a message if no peaks are found.

    Raises:
        KeyError: If the expected dataset structure is not found in the HDF5 file.
    """

    expts_with_peaks = {}
    with h5py.File(file_path, "r") as fp:
        for expt_id in fp["JasonDocument/NMR/NMRData"]:
            try:
                peakList = fp[f"JasonDocument/NMR/NMRData/{expt_id}/Peaks/PeakList"]
                print(len(peakList), "peaks found in", expt_id)
                if len(peakList) > 0:
                    expts_with_peaks[expt_id] = {}
                    expts_with_peaks[expt_id]["spec_info"] = {}
                    expts_with_peaks[expt_id]["peaks"] = {}
                    expts_with_peaks[expt_id]["integrals"] = {}
                    expts_with_peaks[expt_id]["multiplets"] = {}
            except KeyError:
                print(f"No peaks found in {expt_id}")
                continue
    return expts_with_peaks


def get_spec_info(fn: Path, expt_id: str) -> dict:
    """
    Extracts and returns NMR experiment specification information from a JEOL HDF5 file.

    Args:
        fn (Path): Path to the HDF5 file containing NMR experiment data.
        expt_id (str): Identifier for the specific experiment within the file.

    Returns:
        dict: A dictionary containing extracted specification information, including:
            - pulsesequence (str): Pulse sequence used in the experiment.
            - experimenttype (str): Type of experiment (e.g., 'H1_1D', 'C13_1D', 'COSY').
            - datafilename (str): Original data file name.
            - solvent (str): Solvent used in the experiment.
            - expt_fn (str): Filename of the experiment data (without directory).
            - specfrequency (list): List of spectrometer frequencies (positive values).
            - type (str): Dimensionality of the experiment (e.g., '1D', '2D').
            - temperature (float or None): Temperature at which the experiment was performed.
            - nucleus (str or list): Nucleus or nuclei involved in the experiment.

    Notes:
        - Handles decoding of string attributes where necessary.
        - Adjusts experiment type for certain cases based on nucleus and initial type.
        - Handles missing keys gracefully by printing an error message.
    """

    with h5py.File(fn, "r") as fp:

        specInfo = {}
        try:
            spec_info = fp[f"JasonDocument/NMR/NMRData/{expt_id}/SpecInfo"]
            nucleides = fp[f"JasonDocument/NMR/NMRData/{expt_id}/SpecInfo/Nucleides"]

            # Decode if attribute has decode function
            for key, attr in zip(
                ["pulsesequence", "experimenttype", "datafilename", "solvent"],
                [
                    "PulseProgram",
                    "ExperimentType.str",
                    "OrigFilename.filepath.str",
                    "Solvent",
                ],
            ):
                value = spec_info.attrs.get(attr, "")
                if hasattr(value, "decode"):
                    value = value.decode()
                specInfo[key] = value

            # add the filename alone not all the directory
            specInfo["expt_fn"] = Path(specInfo["datafilename"]).name

            specInfo["specfrequency"] = [
                sf
                for sf in spec_info.attrs.get("SpectrometerFrequencies", [])
                if sf > 0
            ]
            specInfo["type"] = (
                str(len(specInfo["specfrequency"])) + "D"
            )  # Number of frequencies
            specInfo["temperature"] = spec_info.attrs.get("Temperature", None)

            Isotopes = [nucleides[n].attrs["Isotope"] for n in list(nucleides)]
            Names = [nucleides[n].attrs["Name"].decode() for n in list(nucleides)]
            nuclei = [f"{i}{n}" for n, i in zip(Names, Isotopes) if i > 0]
            if len(nuclei) == 1:
                specInfo["nucleus"] = nuclei[0]
            else:
                specInfo["nucleus"] = nuclei

            # It seems JEOL cannot figure out H1_1D and C13_1D dataset so we need to set the experiment type manually
            #  based on the dataset being 1D and the nucleus and the experiement type intiallly being set to 'Conventional pulse acquire'
            if specInfo["experimenttype"] == "Conventional pulse acquire":
                if "1H" == specInfo["nucleus"]:
                    specInfo["experimenttype"] = "H1_1D"
                elif "13C" == specInfo["nucleus"]:
                    specInfo["experimenttype"] = "C13_1D"

            # jeol, COSY uses "generic COSY" so we need to change it to "COSY"
            if specInfo["experimenttype"] == "generic COSY":
                specInfo["experimenttype"] = "COSY"

            # datasets[expt_id]["specInfo"] = specInfo
        except KeyError:
            print(f"KeyError: SpecInfo not found for experiment {expt_id}")

    return specInfo


def get_peaklist(filename: Path, expt_id: str) -> list:
    """
    Retrieves the peak list from an HDF5 file for a specified NMR experiment.

    Args:
        filename (Path): Path to the HDF5 file containing NMR data.
        expt_id (str): Identifier for the NMR experiment.

    Returns:
        list: A list of peaks from the specified experiment. Returns an empty list if the peak list is not found.
    """
    try:
        with h5py.File(filename, "r") as fp:
            return list(fp[f"JasonDocument/NMR/NMRData/{expt_id}/Peaks/PeakList"])
    except KeyError:
        return []


def get_peakinfo(filename: Path, expt_id: str, pk_id: str) -> dict:
    """
    Extracts peak information from an HDF5 file for a given experiment and peak ID.

    Parameters:
        filename (Path): Path to the HDF5 file containing NMR data.
        expt_id (str): Experiment ID within the HDF5 file.
        pk_id (str): Peak ID within the experiment's peak list.

    Returns:
        dict: A dictionary containing peak information with keys:
            - "intensity": The peak height (float or None).
            - "delta2": The first position value (float or None).
            - "delta1": The second position value (float or None).
            - "annotation": An empty string (reserved for future use).
            - "type": An integer (default 0).
        Returns an empty dictionary if the specified peak is not found.
    """
    peak_vals = {}
    try:
        with h5py.File(filename, "r") as fp:
            peak_group = fp[
                f"JasonDocument/NMR/NMRData/{expt_id}/Peaks/PeakList/{pk_id}"
            ]
            peak_vals["intensity"] = peak_group.attrs.get("Height", None)
            peak_vals["delta2"] = peak_group.attrs.get("Pos", None)[0]
            peak_vals["delta1"] = peak_group.attrs.get("Pos", None)[1]
            peak_vals["annotation"] = ""
            peak_vals["type"] = 0

        return peak_vals
    except KeyError:
        return {}


def get_integrallist(filename: Path, expt_id: str) -> list:
    """
    Retrieves the list of integrals for a specified experiment from an HDF5 file.

    Args:
        filename (Path): Path to the HDF5 file containing NMR data.
        expt_id (str): Experiment ID used to locate the relevant dataset within the file.

    Returns:
        list: A list of integrals from the specified experiment. Returns an empty list if the dataset is not found.
    """
    try:
        with h5py.File(filename, "r") as fp:
            return list(
                fp[
                    f"JasonDocument/NMR/NMRData/{expt_id}/Multiplets_Integrals/MultipletList"
                ]
            )
    except KeyError:
        return []


def get_integralinfo(filename: Path, expt_id: str, int_id: str) -> dict:
    """
    Extracts integral information from a JEOL NMR HDF5 file for a specified experiment and integral ID.

    Parameters
    ----------
    filename : Path
        Path to the HDF5 file containing NMR data.
    expt_id : str
        Experiment identifier within the HDF5 file.
    int_id : str
        Integral identifier within the experiment.

    Returns
    -------
    dict
        Dictionary containing the following keys:
            - "intensity": The integral value (float or None).
            - "rangeMin1": Minimum value of the first spectrum range (float or None).
            - "rangeMax1": Maximum value of the first spectrum range (float or None).
            - "rangeMin2": Minimum value of the second spectrum range (float or None).
            - "rangeMax2": Maximum value of the second spectrum range (float or None).
            - "delta1": Midpoint of the first spectrum range (float or None).
            - "delta2": Midpoint of the second spectrum range (float or None).
            - "annotation": Annotation string (always empty).
            - "type": Type indicator (always 0).
        Returns an empty dictionary if the specified group is not found.
    """
    integral_vals = {}
    try:
        with h5py.File(filename, "r") as fp:
            integral_group = fp[
                f"JasonDocument/NMR/NMRData/{expt_id}/Multiplets_Integrals/MultipletList/{int_id}"
            ]
            integral_vals["intensity"] = integral_group.attrs.get("Value", None)
            integral_vals["rangeMin1"] = integral_group.attrs.get(
                "SpectrumRange[0]", None
            )[0]
            integral_vals["rangeMax1"] = integral_group.attrs.get(
                "SpectrumRange[0]", None
            )[1]
            integral_vals["rangeMin2"] = integral_group.attrs.get(
                "SpectrumRange[1]", None
            )[0]
            integral_vals["rangeMax2"] = integral_group.attrs.get(
                "SpectrumRange[1]", None
            )[1]
            integral_vals["delta1"] = (
                (integral_vals["rangeMin1"] + integral_vals["rangeMax1"]) / 2
                if integral_vals["rangeMin1"] is not None
                and integral_vals["rangeMax1"] is not None
                else None
            )
            integral_vals["delta2"] = (
                (integral_vals["rangeMin2"] + integral_vals["rangeMax2"]) / 2
                if integral_vals["rangeMin2"] is not None
                and integral_vals["rangeMax2"] is not None
                else None
            )
            integral_vals["annotation"] = ""
            integral_vals["type"] = 0
            return integral_vals
    except KeyError:
        return {}


def createWorkingDirectory(fn: Path) -> dict:
    """
    Creates a working directory structure based on the parent directory of the given file path.
    Args:
        fn (Path): A pathlib.Path object representing the file path.
    Returns:
        dict: A dictionary containing information about the working directory, including its type,
              count, and the directory path with forward slashes.
    """
    # Create the working directory structure

    print("Creating working directory for:", fn)

    workingDir = fn.parent
    print("Working directory is:", workingDir)
    workingDir = str(workingDir).replace("\\", "/")
    print("Working directory with forward slashes:", workingDir)

    working_dir = {
        "datatype": "workingDirectory",
        "count": 1,
        "data": {"0": workingDir},
    }

    return working_dir


def createWorkingFilename(fn: Path) -> dict:
    """
    Generates a dictionary representing a working filename from a given Path object.
    Args:
        fn (Path): The Path object representing the file.
    Returns:
        dict: A dictionary with the following structure:
            {
                "data": {"0": <filename_without_extension>}
            where <filename_without_extension> is the stem of the provided file name.
    """

    working_filename = {
        "datatype": "workingFilename",
        "count": 1,
        "data": {"0": fn.name.split(".")[0]},
    }

    return working_filename


# List of chemical element symbols indexed by atomic number (starting from 0)
symbolElements = (
    "",  # 0: No element
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
    "Rg",
    "Cn",
    "Nh",
    "Fl",
    "Mc",
    "Lv",
    "Ts",
    "Og",
)


def create_carbonatomsinfo(atoms_df: pd.DataFrame) -> dict:
    """
    Generates a dictionary containing information about carbon atoms from a given DataFrame.
    The function iterates over the rows of the input DataFrame, identifies carbon atoms (where the element symbol is 'C'),
    and collects relevant information for each carbon atom. The resulting dictionary includes metadata and a mapping of
    atom indices to their respective information.
    Args:
        atoms_df (pd.DataFrame): A pandas DataFrame containing atom data. Each row should represent an atom and include
            at least the columns 'El' (element index) and 'nH' (number of protons/hydrogens).
    Returns:
        dict: A dictionary with the following structure:
            {
                "count": <number of carbon atoms>,
                "data": {
                    <atom_idx>: {
                        "atom_idx": <index in DataFrame>,
                        "id": <index in DataFrame>,
                        "atomNumber": <index + 1>,
                        "symbol": "C",
                        "numProtons": <number of protons/hydrogens>
                    },
                    ...
    """

    carbon_atoms_info = {
        "datatype": "carbonAtomsInfo",
        "count": atoms_df.shape[0],
        "data": {},
    }

    for idx, row in atoms_df.iterrows():
        symbol = symbolElements[row.get("El", 0)]
        if symbol != "C":
            carbon_atoms_info["count"] -= 1
            continue
        atom_info = {
            "atom_idx": idx,
            "id": idx,
            "atomNumber": idx + 1,
            "symbol": symbol,
            "numProtons": row.get("nH", None),
        }
        carbon_atoms_info["data"][str(idx)] = atom_info

    return carbon_atoms_info


def create_allatomsinfo(atoms_df: pd.DataFrame) -> dict:
    """
    Generates a dictionary containing information about all atoms from a given DataFrame.

    Args:
        atoms_df (pd.DataFrame): A pandas DataFrame containing atom data. Each row should represent an atom and include
            at least the columns 'El' (element index) and 'nH' (number of protons/hydrogens).

    Returns:
        dict: A dictionary with the following structure:
            {
                "datatype": "allAtomsInfo",
                "count": <number of atoms>,
                "data": {
                    <atom_idx>: {
                        "atom_idx": <index in DataFrame>,
                        "id": <index in DataFrame>,
                        "atomNumber": <index + 1>,
                        "symbol": <element symbol>,
                        "numProtons": <number of protons/hydrogens>
                    },
                    ...
                }
            }
    """
    all_atoms_info = {
        "datatype": "allAtomsInfo",
        "count": atoms_df.shape[0],
        "data": {},
    }

    for idx, row in atoms_df.iterrows():
        symbol = symbolElements[row.get("El", 0)]
        atom_info = {
            "atom_idx": idx,
            "id": idx,
            "atomNumber": idx + 1,
            "symbol": symbol,
            "numProtons": row.get("nH", None),
        }
        all_atoms_info["data"][str(idx)] = atom_info

    return all_atoms_info


def readJEOLmolecule(fn: Path) -> pd.DataFrame:
    """
    Reads atom information from a JEOL HDF5 molecule file and returns it as a pandas DataFrame.

    Args:
        fn (Path): Path to the JEOL HDF5 file.

    Returns:
        pd.DataFrame: DataFrame containing atom information with columns for atom attributes and a 'Z' axis (set to 0.0).
                      The DataFrame is sorted by 'atom_idx' and indexed from 0.

    Notes:
        - If the required HDF5 keys are missing, a KeyError is caught and printed.
        - The function assumes a specific structure for the JEOL HDF5 file.
    """
    with h5py.File(fn, "r") as fp:
        try:
            mol_data = fp["JasonDocument/Molecules/Molecules/0/Atoms"]
            print(mol_data.keys())
            atom_list = []
            for atom_idx in mol_data.keys():
                atom_info = {
                    key: mol_data[atom_idx].attrs.get(key, None)
                    for key in mol_data[atom_idx].attrs.keys()
                }
                atom_info["atom_idx"] = int(atom_idx)
                atom_list.append(atom_info)
            atoms_df = pd.DataFrame(atom_list)
            atoms_df["Z"] = 0.0
        except KeyError as e:
            print(f"KeyError: {e}")

    atoms_df = atoms_df.sort_values(by="atom_idx").reset_index(drop=True)
    return atoms_df


def createMNOVAcalcMethod() -> dict:
    """
    Creates a dictionary representing the MNOVA calculation method.

    Returns:
        dict: A dictionary with the following structure:
            {
                "datatype": "MNOVAcalcMethod",
                "count": 1,
                "data": {"0": "JEOL Predict"}
            }
    """
    mNOVA_calc_method = {
        "datatype": "MNOVAcalcMethod",
        "count": 1,
        "data": {"0": "JEOL Predict"},
    }
    return mNOVA_calc_method


def createCarbonCalcPositionsMethod():
    def createCarbonCalcPositionsMethod() -> dict:
        """
        Creates and returns a dictionary representing the carbon calculation positions method structure.

        Returns:
            dict: A dictionary with the following structure:
                {
        """

    carbon_calc_positions_method = {
        "datatype": "carbonCalcPositionsMethod",
        "count": 1,
        "data": {"0": "Calculated Positions"},
    }
    return carbon_calc_positions_method


def dataframe_to_rdkit_molecule(
    df: pd.DataFrame,
) -> Chem.Mol:
    """
    Convert a pandas DataFrame containing molecular structure data to an RDKit molecule.

    Args:
        df (pd.DataFrame): DataFrame with columns ['El', 'NB.Conn', 'NB.Num', 'X', 'Y', 'nH', 'Z'].
            - 'El': Atomic number of each atom.
            - 'NB.Conn': List of bond orders for each connection.
            - 'NB.Num': List of connected atom indices.
            - 'X', 'Y', 'Z': Cartesian coordinates for each atom.
            - 'nH': Number of hydrogens (optional).

    Returns:
        Chem.Mol: RDKit molecule object constructed from the DataFrame, or None if creation failed.

    Notes:
        - Bonds are added based on 'NB.Num' and 'NB.Conn' columns.
        - Atom coordinates are set if 'X' and 'Y' columns are present.
        - Attempts to sanitize the molecule; returns unsanitized molecule if sanitization fails.
    """

    try:
        # Create an editable molecule object
        mol = Chem.EditableMol(Chem.Mol())

        # Reset index to ensure proper atom numbering
        mol_df = df.copy().reset_index(drop=True)

        # Handle missing Z coordinates
        if "Z" not in mol_df.columns:
            mol_df["Z"] = 0.0
        mol_df["Z"] = mol_df["Z"].fillna(0.0)

        # Step 1: Add all atoms
        atom_indices: dict[int, int] = {}
        for idx, row in mol_df.iterrows():
            element_num = int(row["El"])
            atom = Chem.Atom(element_num)
            atom_idx = mol.AddAtom(atom)
            atom_indices[idx] = atom_idx

        # Step 2: Add bonds based on connectivity
        bonds_added: set[tuple[int, int]] = set()
        bond_count = 0

        for idx, row in mol_df.iterrows():
            rdkit_atom_idx = atom_indices[idx]
            nb_num = row.get("NB.Num")
            nb_conn = row.get("NB.Conn")

            nb_num_valid = nb_num is not None and not (
                isinstance(nb_num, float) and pd.isna(nb_num)
            )
            nb_conn_valid = nb_conn is not None and not (
                isinstance(nb_conn, float) and pd.isna(nb_conn)
            )

            if nb_num_valid and nb_conn_valid:
                try:
                    if isinstance(nb_num, np.ndarray):
                        nb_num = nb_num.tolist()
                    elif not isinstance(nb_num, list):
                        nb_num = [nb_num]

                    if isinstance(nb_conn, np.ndarray):
                        nb_conn = nb_conn.tolist()
                    elif not isinstance(nb_conn, list):
                        nb_conn = [nb_conn]

                    for i, connected_atom in enumerate(nb_num):
                        if (
                            connected_atom is not None
                            and not (
                                isinstance(connected_atom, float)
                                and pd.isna(connected_atom)
                            )
                            and connected_atom != ""
                        ):
                            try:
                                connected_df_idx = int(float(connected_atom))
                                if connected_df_idx == idx:
                                    continue
                                if connected_df_idx not in atom_indices:
                                    print(
                                        f"   Warning: Connected atom {connected_df_idx} not found in dataframe"
                                    )
                                    continue
                                connected_rdkit_idx = atom_indices[connected_df_idx]
                                bond_pair = tuple(
                                    sorted([rdkit_atom_idx, connected_rdkit_idx])
                                )
                                if bond_pair not in bonds_added:
                                    if i < len(nb_conn):
                                        bond_order = nb_conn[i]
                                    else:
                                        bond_order = 1

                                    if bond_order == 1:
                                        rdkit_bond_type = Chem.BondType.SINGLE
                                    elif bond_order == 2:
                                        rdkit_bond_type = Chem.BondType.DOUBLE
                                    elif bond_order == 3:
                                        rdkit_bond_type = Chem.BondType.TRIPLE
                                    elif bond_order in (513, 770, 769):
                                        rdkit_bond_type = Chem.BondType.AROMATIC
                                    elif bond_order == 514:
                                        rdkit_bond_type = Chem.BondType.DOUBLE
                                    else:
                                        rdkit_bond_type = Chem.BondType.SINGLE

                                    mol.AddBond(
                                        rdkit_atom_idx,
                                        connected_rdkit_idx,
                                        rdkit_bond_type,
                                    )
                                    bonds_added.add(bond_pair)
                                    bond_count += 1

                            except (ValueError, TypeError) as e:
                                print(
                                    f"   Warning: Could not parse connection {connected_atom} for atom {idx}: {e}"
                                )
                                continue

                except Exception as e:
                    print(
                        f"   Warning: Error processing connections for atom {idx}: {e}"
                    )
                    continue

        # Step 3: Convert to final molecule and add 2D coordinates
        final_mol = mol.GetMol()

        if all(col in mol_df.columns for col in ["X", "Y"]):
            conf = Chem.Conformer(len(mol_df))
            conf.Set3D(False)
            for idx, row in mol_df.iterrows():
                rdkit_idx = atom_indices[idx]
                x, y = float(row["X"]), float(row["Y"])
                conf.SetAtomPosition(rdkit_idx, (x, y, 0.0))
            final_mol.AddConformer(conf)
        else:
            print(
                "X and/or Y coordinates not found - molecule will have no coordinates"
            )

        # Step 4: Sanitize the molecule
        try:
            Chem.SanitizeMol(final_mol)
            return final_mol

        except Exception as sanitize_error:
            print(f"Sanitization failed: {sanitize_error}")
            print("Trying with less strict sanitization...")

            try:
                sanitize_ops = (
                    Chem.SanitizeFlags.SANITIZE_FINDRADICALS
                    | Chem.SanitizeFlags.SANITIZE_KEKULIZE
                    | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY
                    | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION
                    | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION
                    | Chem.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY
                )
                Chem.SanitizeMol(final_mol, sanitizeOps=sanitize_ops)
                print("Molecule sanitized with custom operations!")
                return final_mol

            except Exception as e2:
                print(f"Custom sanitization also failed: {e2}")
                print("Returning unsanitized molecule...")
                return final_mol

    except Exception as e:
        print(f"Error creating molecule: {e}")
        return None


def fix_stereo_info(molstr):
    def fix_stereo_info(molstr: str) -> str:
        """
        Fixes stereo information in a molecular string representation.

        This function processes a multi-line string representing a molecule, where each line contains
        atom connectivity information. For lines with exactly four columns, if the last column (stereo info)
        is not '0', it replaces it with '0' (or '00' for two-character stereo info) to standardize the stereo information.

        Args:
            molstr (str): The molecular string to be processed.

        Returns:
            str: The modified molecular string with fixed stereo information.
        """

    lines = molstr.splitlines()
    for i, line in enumerate(lines):
        words = line.split()
        # Check if the line has 4 columns then we have reached the connectivity group
        if len(words) == 4:
            # If the last column is not '0', we have stereo info to fix
            if words[-1] != "0":
                print(f"Fixing stereo info in line {i}: {line}")
                # Fix the stereo info by replacing it with '0'
                if len(words[-1]) == 1:
                    line = line[:-1] + "0"
                elif len(words[-1]) == 2:
                    line = line[:-2] + "00"
                else:
                    print("Unexpected length:", len(words[-1]))
                # update the line in the original list
                lines[i] = line

    return "\n".join(lines)


def create_molfile(molstr: str) -> dict:
    """
    Creates a dictionary representing a molfile with the given molecular string.
    Args:
        molstr (str): A string containing the molecular data in molfile format.
    Returns:
        dict: A dictionary with keys 'datatype', 'count', and 'data', where 'data' contains the molfile string.
    """

    molfile = {"datatype": "molfile", "count": 1, "data": {"0": molstr}}
    return molfile


def create_smiles(smiles_str: str) -> dict:
    """
    Creates a dictionary representing a SMILES string.
    Args:
        smiles_str (str): A string containing the SMILES representation of a molecule.
    Returns:
        dict: A dictionary with keys 'datatype', 'count', and 'data', where 'data' contains the SMILES string.
    """
    smiles = {"datatype": "smiles", "count": 1, "data": {"0": smiles_str}}
    return smiles


def get_chosenSpectra(specAssignments: list) -> dict:
    """
    Generates a dictionary summarizing selected spectra assignments.

    Args:
        specAssignments (list): A list of dictionaries, each containing
            'experiment_name' and 'experiment_type' keys representing individual spectra assignments.

    Returns:
        dict: A dictionary with the following structure:
            {
                "count": <number of assignments>,
                "data": {
                    "<index>": "<experiment_name> <experiment_type>",
                    ...
        The 'data' field maps string indices to a concatenation of experiment name and type.
    """

    chosenSpectra = {
        "datatype": "chosenSpectra",
        "count": len(specAssignments),
        "data": {},
    }
    for i, data in enumerate(specAssignments):
        chosenSpectra["data"][
            str(i)
        ] = f"{data['experiment_name']} {data['experiment_type']}"

    return chosenSpectra


def get_spectraWithPeaks(specAssignments: list) -> dict:
    """
    Generates a dictionary summarizing spectra assignments with their experiment names.

    Args:
        specAssignments (list): A list of dictionaries, each containing spectral assignment data.
            Each dictionary is expected to have an "experiment_name" key.

    Returns:
        dict: A dictionary with the following structure:
            {
                "count": <number of assignments>,
                "data": {
                    "<index>": <experiment_name>,
                    ...
    """

    spectraWithPeaks = {
        "datatype": "spectraWithPeaks",
        "count": len(specAssignments),
        "data": {},
    }
    for i, data in enumerate(specAssignments):
        spectraWithPeaks["data"][str(i)] = data["experiment_name"]

    return spectraWithPeaks


class jeolData:
    """
    A class for handling and analyzing JEOL NMR data files, extracting molecular and spectral information,
    and preparing data for further processing and export.
    Attributes:
        analysis_cancelled (bool): Indicates if the analysis was cancelled.
        spectra_assignments (dict): Stores assignments of spectra to experiment types.
        file_path (Path): Path to the JEOL data file.
        spectra_used_in_analysis (dict): Spectra used in the analysis.
        chosenSpectra (dict): Spectra chosen for analysis.
        spectraWithPeaks (dict): Spectra containing peak information.
        workingDirectory (str): Working directory for analysis.
        workingFilename (str): Working filename for analysis.
        datasets_with_peaks (dict): Datasets containing peak information.
        atoms_df (DataFrame): DataFrame containing atom information from the molecule.
        allatomsInfo (dict): Information about all atoms in the molecule.
        carbonatomsInfo (dict): Information about carbon atoms in the molecule.
        rdkit_mol (Mol): RDKit molecule object.
        molstr (str): Molfile string with stereo information fixed.
        molfile (str): Final molfile string.
        smiles_str (str): SMILES string representation of the molecule.
        smiles (str): Final SMILES string.
        hostname (str): Hostname of the machine running the analysis.
        c13_predictions_df (DataFrame): DataFrame of C13 predictions.
        c13predictions (dict): Processed C13 prediction data.
        MNOVAcalcMethod (str): Calculation method for MNOVA.
        carbonCalcPositionsMethod (str): Calculation method for carbon positions.
        simulated_annealing (Any): Simulated annealing processing option.
        ml_consent (Any): Machine learning consent option.
    Methods:
        updateDatawithExptNames():
            Prints experiment names and associated metadata for datasets with peaks.
        createJsonDict():
            Creates a dictionary suitable for JSON export containing all relevant molecular and spectral data.
        append_specdata():
            Appends selected spectral data to the export dictionary, organizing by experiment type.
        choosePeakPickedSpectaforSimpleNMR():
            Presents a dialog for the user to select and assign experiment types to available spectra,
            and updates internal data structures accordingly.
    """

    def __init__(self, file_path: Path):
        """
        Initialize the analysis object for JEOL NMR data.

        Parameters
        ----------
        file_path : Path
            Path to the JEOL NMR data file.

        Attributes
        ----------
        analysis_cancelled : bool
            Flag indicating if the analysis was cancelled.
        spectra_assignments : dict
            Dictionary to store spectra assignments.
        file_path : Path
            Path to the input data file.
        spectra_used_in_analysis : dict
            Spectra used in the analysis.
        chosenSpectra : dict
            Chosen spectra for analysis.
        spectraWithPeaks : dict
            Spectra containing peak information.
        workingDirectory : str
            Working directory for analysis files.
        workingFilename : str
            Working filename for analysis files.
        datasets_with_peaks : dict
            Datasets containing peak information.
        atoms_df : pandas.DataFrame
            DataFrame containing atom information from JEOL molecule file.
        allatomsInfo : dict
            Information about all atoms in the molecule.
        carbonatomsInfo : dict
            Information about carbon atoms in the molecule.
        rdkit_mol : rdkit.Chem.Mol
            RDKit molecule object.
        molstr : str
            Molfile string representation of the molecule.
        molfile : str
            Final molfile string after processing.
        smiles_str : str
            SMILES string representation of the molecule.
        smiles : str
            Final SMILES string after processing.
        hostname : str
            Hostname of the machine running the analysis.
        c13_predictions_df : pandas.DataFrame
            DataFrame containing C13 prediction data.
        c13predictions : dict
            Processed C13 prediction information.
        MNOVAcalcMethod : str
            Calculation method for MNOVA.
        carbonCalcPositionsMethod : str
            Calculation method for carbon positions.

        Notes
        -----
        This initializer reads and processes JEOL NMR data, extracts molecular and spectral information,
        and prepares all necessary attributes for further analysis.
        """

        self.analysis_cancelled = False
        self.spectra_assignments = {}
        self.file_path = file_path

        self.spectra_used_in_analysis = {}
        self.chosenSpectra = {}
        self.spectraWithPeaks = {}

        self.workingDirectory = createWorkingDirectory(file_path)
        self.workingFilename = createWorkingFilename(file_path)

        self.datasets_with_peaks = find_datasets_with_peaks(file_path)

        #  add allatomsinfo
        self.atoms_df = readJEOLmolecule(self.file_path)
        self.allatomsInfo = create_allatomsinfo(self.atoms_df)

        # add carbonatomsinfo
        self.carbonatomsInfo = create_carbonatomsinfo(self.atoms_df)
        print(self.carbonatomsInfo)

        # create rdkit mol object
        self.rdkit_mol = dataframe_to_rdkit_molecule(self.atoms_df)

        # create molfile str
        self.molstr = Chem.MolToMolBlock(self.rdkit_mol)
        self.molstr = fix_stereo_info(self.molstr)
        self.molfile = create_molfile(self.molstr)

        # create smiles string
        self.smiles_str = Chem.MolToSmiles(self.rdkit_mol)

        self.smiles = create_smiles(self.smiles_str)

        # add hostname
        self.hostname = get_host_name()

        # get c13 predictions
        self.c13_predictions_df = extractC13predictionsJEOL(
            self.file_path, self.atoms_df
        )

        self.c13predictions = create_c13predictions(self.c13_predictions_df)

        # add spec_info to data with peaks
        for expt_id, data in self.datasets_with_peaks.items():
            data["spec_info"] = get_spec_info(self.file_path, expt_id)

        # add peakList to data with peaks
        for expt_id, data in self.datasets_with_peaks.items():
            data["peaks"]["data"] = {}
            data["peaks"]["count"] = 0
            data["peaks"]["datatype"] = "peaks"
            for pk_id in get_peaklist(self.file_path, expt_id):
                peak_info = get_peakinfo(self.file_path, expt_id, pk_id)
                data["peaks"]["data"][pk_id] = peak_info
                data["peaks"]["count"] += 1

        # Add integral information to data with peaks
        for expt_id, data in self.datasets_with_peaks.items():
            data["integrals"]["data"] = {}
            data["integrals"]["count"] = 0
            data["integrals"]["normValue"] = 1
            data["integrals"]["datatype"] = "integrals"
            for int_id in get_integrallist(self.file_path, expt_id):
                integral_info = get_integralinfo(self.file_path, expt_id, int_id)
                data["integrals"]["data"][int_id] = integral_info
                data["integrals"]["count"] += 1

        # Add multiplet information to data with peaks
        # for now we will just create empty multiplet info
        for expt_id, data in self.datasets_with_peaks.items():
            data["multiplets"]["data"] = {}
            data["multiplets"]["count"] = 0
            data["multiplets"]["normValue"] = 1
            data["multiplets"]["datatype"] = "multiplets"

        # add workingDirectory
        self.workingDirectory = createWorkingDirectory(self.file_path)

        # add workingFilename
        self.workingFilename = createWorkingFilename(self.file_path)

        # add MNOVAcalcMethod
        self.MNOVAcalcMethod = createMNOVAcalcMethod()

        # add carbonCalcPositionsMethod
        self.carbonCalcPositionsMethod = createCarbonCalcPositionsMethod()

    def createJsonDict(self) -> dict:
        """
        Creates and returns a dictionary containing all relevant data for JSON serialization.

        The dictionary includes molecular information, calculation methods, atom data, predictions,
        spectra information, simulated annealing results, machine learning consent, and additional
        spectral data appended from `append_specdata()`.

        Returns:
            dict: A dictionary with all necessary fields for JSON output.
        """

        json_dict = {
            "smiles": self.smiles,
            "molfile": self.molfile,
            "hostname": self.hostname,
            "workingDirectory": self.workingDirectory,
            "workingFilename": self.workingFilename,
            "MNOVAcalcMethod": self.MNOVAcalcMethod,
            "carbonCalcPositionsMethod": self.carbonCalcPositionsMethod,
            "allAtomsInfo": self.allatomsInfo,
            "carbonAtomsInfo": self.carbonatomsInfo,
            "c13predictions": self.c13predictions,
            "chosenSpectra": self.chosenSpectra,
            "spectraWithPeaks": self.spectraWithPeaks,
            "simulatedAnnealing": {
                "datatype": "simulatedAnnealing",
                "count": 1,
                "data": {"0": self.simulated_annealing},
            },
            "ml_consent": {
                "datatype": "ml_consent",
                "count": 1,
                "data": {"0": self.ml_consent},
            },
        }

        spec_dict = self.append_specdata()

        # merge the two dictionaries
        json_dict.update(spec_dict)

        return json_dict

    def append_specdata(self) -> dict:
        """
        Appends selected spectral data to a dictionary for JSON export.

        This method processes the chosen spectra assignments, matches them to datasets with peaks,
        and organizes the spectral data by experiment type and index. It includes information such as
        spec_info, peaks, integrals, and multiplets for each spectrum.

        Returns:
            dict: A dictionary where keys are spectrum identifiers (e.g., 'HSQC_0', 'COSY_1') and values
                  are dictionaries containing spectral data and metadata.
        """
        expts_type_count: dict[str, int] = {
            "SKIP": 0,
            "H1_1D": 0,
            "C13_1D": 0,
            "Pureshift": 0,
            "DEPT": 0,
            "COSY": 0,
            "NOESY": 0,
            "HSQC": 0,
            "HMBC": 0,
            "DDEPT_CH3_ONLY": 0,
            "HSQC_CLIPCOSY": 0,
        }

        expt_kys: list[str] = []
        expt_types: list[str] = []
        for idx, expt_id in self.chosenSpectra["data"].items():
            # split expt_id into a list based on space
            expt_id_parts = expt_id.split(" ")
            # if length of list greater than 4 then combine the first two items together with " "
            if len(expt_id_parts) > 4:
                expt_id_parts = [" ".join(expt_id_parts[:2])] + expt_id_parts[2:]

            expt_ky = expt_id_parts[-2]
            expt_type = expt_id_parts[-1]

            expt_kys.append(expt_ky)
            expt_types.append(expt_type)

        # get the experiment from data_with_peaks based on expt_ky
        exptNameToKey: dict[str, str] = {}

        for id, data in self.datasets_with_peaks.items():
            exptNameToKey[data["spec_info"]["expt_fn"]] = id
            print(f"Experiment: {data['spec_info']['expt_fn']}, ID: {id}")

        spectra: dict[str, dict] = {}

        for expt_ky, expt_type in zip(expt_kys, expt_types):
            data = None

            data = self.datasets_with_peaks.get(expt_ky, None)
            if data:
                spec_id = f"{expt_type}_{expts_type_count[expt_type]}"
                spectra[spec_id] = {}
                expts_type_count[expt_type] += 1
                spectra[spec_id]["datatype"] = "nmrspectrum"

                # add spec_info
                for k, v in data["spec_info"].items():
                    spectra[spec_id][k] = v
                # add peaks
                spectra[spec_id]["peaks"] = data["peaks"]
                # add integrals
                spectra[spec_id]["integrals"] = data["integrals"]
                # add multiplets
                spectra[spec_id]["multiplets"] = data["multiplets"]
            else:
                print(f"Experiment {expt_ky} not found")

        return spectra

    def choosePeakPickedSpectaforSimpleNMR(self):
        """
        Presents a dialog for the user to select peak-picked spectra for further NMR analysis.

        This method generates a list of experiment names from datasets with peaks, displays them in a dialog,
        and allows the user to assign experiment types or skip certain spectra. Based on the user's selections,
        it updates the spectra assignments, chosen spectra, and spectra with peaks. It also retrieves processing
        options such as simulated annealing and machine learning consent from the dialog.

        If the dialog is cancelled by the user, the analysis is marked as cancelled.

        Returns:
            None
        """

        # Generate random experiment names (simulating found experiments)
        experiment_names = []
        chosen_types = {}

        # create experiment_names info from self.datasets_with_peaks information

        for id, data in self.datasets_with_peaks.items():

            expt_type = data["spec_info"]["experimenttype"]
            expt_fn = Path(data["spec_info"]["datafilename"]).name
            expt_nucleus = data["spec_info"]["nucleus"].__str__()
            # remove quotes from expt_nucleus
            expt_nucleus = expt_nucleus.replace('"', "")
            expt_nucleus = expt_nucleus.replace("'", "")

            expt_name = (
                f'{expt_nucleus} {data["spec_info"]["pulsesequence"]} {expt_fn} {id}'
            )
            experiment_names.append(expt_name)
            chosen_types[expt_name] = expt_type

        # Create and show the dialog
        dialog = NMRExperimentDialog(experiment_names, chosen_types=chosen_types)

        # Execute dialog and get result
        result = dialog.exec_()

        if result == QDialog.Accepted:
            self.spectra_assignments = dialog.get_experiment_assignments()

            for assignment in self.spectra_assignments:
                # remove from list if assigned type is SKIP
                if assignment["experiment_type"] == "SKIP":
                    self.spectra_assignments.remove(assignment)

            # create chosenSpectra data
            self.chosenSpectra = get_chosenSpectra(self.spectra_assignments)

            # create specrtraWithPeaks data
            self.spectraWithPeaks = get_spectraWithPeaks(self.spectra_assignments)

            self.simulated_annealing, self.ml_consent = dialog.get_processing_options()

        else:
            print("\nDialog was cancelled by user.")
            self.analysis_cancelled = True


class SubmissionOutcome(Enum):
    """Classifies every distinct way /simpleMNOVA can respond.

    See routes.py::simpleMNOVA_display_molecule for the source of truth:
      - SUCCESS: status 200, HTML body -> the rendered D3 result page.
      - DIAGNOSTIC_HTML: non-200 status, but the body is still a full HTML
        document (e.g. a "misassignment of molecule" diagnostic report
        with comparison tables, returned by some pipeline failures with
        status 400 so the user can see where assignment got stuck). This
        is content meant to be viewed, not a plain error string - status
        code alone doesn't distinguish it from ERROR, the body shape does.
      - UNREGISTERED / REGISTRATION_EXPIRED: status 200, but a JSON body
        (Content-Type: application/json) with {"status": ...}. Status code
        alone does NOT distinguish this from SUCCESS - Content-Type does.
      - ERROR: non-200 status with a body that is NOT a full HTML document
        (short plain-text validation errors, or a 200 JSON body that
        doesn't match the known unregistered/registration_expired shape).
      - NETWORK_ERROR: requests.RequestException - no HTTP response at all.
    """

    SUCCESS = "success"
    DIAGNOSTIC_HTML = "diagnostic_html"
    UNREGISTERED = "unregistered"
    REGISTRATION_EXPIRED = "registration_expired"
    ERROR = "error"
    NETWORK_ERROR = "network_error"


@dataclass
class SubmissionResult:
    """Outcome of submit_to_server(), replacing the old bare bool return.

    A bool can't distinguish "unregistered" from "pipeline error" from
    "network unreachable" - each needs different follow-up behaviour in
    __main__ (open the viewer vs. open a registration page vs. just stop),
    so callers must branch on `outcome` rather than truthiness.
    """

    outcome: SubmissionOutcome
    html_path: Optional[Path] = None
    message: Optional[str] = None
    status_code: Optional[int] = None
    registration_url: Optional[str] = None
    log_path: Optional[Path] = None


def _log_submission_error(
    json_data: Dict, status_code, body: str, content_type: str = ""
) -> Optional[Path]:
    """Append raw error details to a log file next to the HTML output.

    Kept separate from the JSON input file so it can accumulate multiple
    submission attempts. Returns the log file path on success, or None if
    logging itself failed (never raises - logging must not mask the
    original error).
    """
    try:
        fn_str = json_data.get("workingDirectory", {}).get("data", {}).get("0", ".")
        log_dir = Path(fn_str, "html")
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / "simpleNMR_submission_errors.log"

        timestamp = datetime.now().isoformat(timespec="seconds")
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(f"\n{'=' * 70}\n")
            f.write(f"Timestamp: {timestamp}\n")
            f.write(f"Status code: {status_code}\n")
            f.write(f"Content-Type: {content_type}\n")
            f.write(f"Body:\n{body}\n")

        return log_path
    except Exception as e:
        print(f"Failed to write submission error log: {e}")
        return None


def _looks_like_html(text: str) -> bool:
    """Sniff whether a response body is a full HTML document.

    Used to distinguish a real page (success result, or a diagnostic
    report like the "misassignment of molecule" tables page, which can
    come back with a non-200 status) from a short plain-text error
    message like "Invalid JSON: ...". Deliberately simple: checks the
    first ~500 characters (skipping any leading HTML comment) for a
    doctype or <html> tag, rather than trying to fully parse the body.
    """
    if not text:
        return False
    head = text.lstrip()[:500].lower()
    return "<!doctype html" in head or "<html" in head


def _write_html_result(json_data: Dict, response_text: str) -> Path:
    """Write a response body to <workingDirectory>/html/<workingFilename>.html.

    Shared by SUCCESS and DIAGNOSTIC_HTML - both need the same file
    written for MainWindow to open, they differ only in the status-bar
    note MainWindow shows once it's open.
    """
    workingFilename = json_data["workingFilename"]["data"].get(
        "0", "nmr_analysis_result"
    )
    response_text = response_text.replace("dummy_title", workingFilename)

    fn_str = json_data["workingDirectory"]["data"].get("0", ".")
    fn_path = Path(fn_str, "html")
    fn_path.mkdir(parents=True, exist_ok=True)
    fn_path = Path(fn_path, workingFilename + ".html")

    with open(fn_path, "w", encoding="utf-8") as f:
        f.write(response_text)

    return fn_path


def submit_to_server(json_data: Dict, simpleNMR_address: str) -> SubmissionResult:
    """
    Submit the JSON data to the processing server with progress dialog.

    Classifies the response per routes.py's actual contract for
    /simpleMNOVA (see SubmissionOutcome docstring) rather than assuming
    status 200 always means "here is the HTML to display" - it can also
    mean an unregistered/expired-registration JSON response.

    Args:
        json_data: The converted JSON data

    Returns:
        SubmissionResult describing what happened and what (if anything)
        the caller should do next.
    """

    # Create progress dialog
    progress = QProgressDialog(
        "Submitting data to simpleNMR server...", "Cancel", 0, 400
    )
    progress.setWindowModality(Qt.WindowModal)
    progress.setMinimumDuration(0)
    progress.setCancelButton(
        None
    )  # Remove cancel button since we can't easily cancel the request
    progress.show()

    # Variables to store result, populated by make_request() in the
    # background thread and read back out once it finishes.
    result = {
        "outcome": None,
        "html_path": None,
        "message": None,
        "status_code": None,
        "registration_url": None,
        "log_path": None,
        "finished": False,
    }

    def make_request():
        try:
            print("Submitting data to simpleNMR server...")

            response = requests.post(
                simpleNMR_address,
                headers={"Content-Type": "application/json"},
                json=json_data,
                timeout=100,
            )

            content_type = response.headers.get("Content-Type", "")
            is_json = "application/json" in content_type.lower()

            if is_json:
                # Could be the unregistered / registration_expired shape,
                # or (defensively) some other JSON we don't recognise.
                try:
                    response_data = response.json()
                except ValueError:
                    response_data = None

                status_field = (
                    response_data.get("status")
                    if isinstance(response_data, dict)
                    else None
                )

                if status_field == "unregistered":
                    result["outcome"] = SubmissionOutcome.UNREGISTERED
                    result["registration_url"] = response_data.get("registration_url")
                    result["message"] = "This machine is not registered."
                    result["status_code"] = response.status_code
                elif status_field == "registration_expired":
                    result["outcome"] = SubmissionOutcome.REGISTRATION_EXPIRED
                    result["registration_url"] = response_data.get("registration_url")
                    result["message"] = "Registration for this machine has expired."
                    result["status_code"] = response.status_code
                else:
                    # Unexpected JSON shape - don't guess, log and surface it.
                    error_msg = f"Unexpected JSON response: {response.text}"
                    print(error_msg)
                    result["outcome"] = SubmissionOutcome.ERROR
                    result["status_code"] = response.status_code
                    result["message"] = error_msg
                    result["log_path"] = _log_submission_error(
                        json_data, response.status_code, response.text, content_type
                    )

            elif _looks_like_html(response.text):
                # A real HTML document, whether or not status is 200.
                # Some pipeline failures deliberately return a full
                # diagnostic report (e.g. a "misassignment of molecule"
                # page with comparison tables) with a non-200 status, so
                # the user can see where assignment got stuck - that's
                # content to display, not a plain error string.
                fn_path = _write_html_result(json_data, response.text)

                if response.status_code == 200:
                    result["outcome"] = SubmissionOutcome.SUCCESS
                else:
                    result["outcome"] = SubmissionOutcome.DIAGNOSTIC_HTML
                    result["status_code"] = response.status_code
                    # Still recorded for the record, even though this
                    # isn't treated as a blocking error - no dialog shown
                    # for this outcome (see below).
                    result["log_path"] = _log_submission_error(
                        json_data, response.status_code, response.text, content_type
                    )

                result["html_path"] = fn_path

            else:
                # Non-200 and not HTML: short plain-text validation errors
                # ("Invalid JSON: ...", "No JSON data received", etc.) or
                # an arbitrary PipelineError message that isn't a report
                # page. Never write/open this as if it were a result.
                error_msg = f"Server error: {response.status_code} - {response.text}"
                print(error_msg)
                result["outcome"] = SubmissionOutcome.ERROR
                result["status_code"] = response.status_code
                result["message"] = response.text
                result["log_path"] = _log_submission_error(
                    json_data, response.status_code, response.text, content_type
                )

        except requests.RequestException as e:
            error_msg = f"Network error: {e}"
            print(error_msg)
            result["outcome"] = SubmissionOutcome.NETWORK_ERROR
            result["message"] = error_msg
        except Exception as e:
            error_msg = f"Error submitting to simpleNMR server: {e}"
            print(error_msg)
            result["outcome"] = SubmissionOutcome.ERROR
            result["message"] = error_msg
            result["log_path"] = _log_submission_error(json_data, None, error_msg)
        finally:
            result["finished"] = True

    # Start request in background thread
    thread = threading.Thread(target=make_request)
    thread.daemon = True
    thread.start()

    # Process events until request is complete
    while not result["finished"]:
        QApplication.processEvents()
        thread.join(0.1)  # Check every 100ms

        # Update progress dialog text periodically to show it's still working
        progress.setValue((progress.value() + 1))

    progress.close()

    outcome = result["outcome"]

    # Show the appropriate dialog for every non-success outcome. SUCCESS
    # is deliberately silent here - the viewer window opening is itself
    # the feedback; __main__ handles opening it.
    if outcome == SubmissionOutcome.NETWORK_ERROR:
        QMessageBox.critical(
            None,
            "Network Error",
            f"Could not reach the simpleNMR server:\n{result['message']}",
        )
    elif outcome == SubmissionOutcome.ERROR:
        log_note = (
            f"\n\nDetails logged to:\n{result['log_path']}"
            if result["log_path"]
            else ""
        )
        QMessageBox.critical(
            None,
            "Submission Error",
            f"Server returned an error"
            f"{' (status ' + str(result['status_code']) + ')' if result['status_code'] is not None else ''}:\n"
            f"{result['message']}{log_note}",
        )
    elif outcome in (
        SubmissionOutcome.UNREGISTERED,
        SubmissionOutcome.REGISTRATION_EXPIRED,
    ):
        title = (
            "Registration Required"
            if outcome == SubmissionOutcome.UNREGISTERED
            else "Registration Expired"
        )
        QMessageBox.warning(None, title, result["message"])
        if result["registration_url"]:
            webbrowser.open(result["registration_url"])

    return SubmissionResult(
        outcome=outcome,
        html_path=result["html_path"],
        message=result["message"],
        status_code=result["status_code"],
        registration_url=result["registration_url"],
        log_path=result["log_path"],
    )


def check_user_registration(address) -> bool:
    """
    Check if the user's machine is registered for the service.

    Returns:
        True if user can proceed, False otherwise
    """
    try:
        # Generate machine ID (MAC address based)
        mac_based_id = hex(uuid.getnode())
        print(f"Machine ID: {mac_based_id}")

        # Prepare request
        json_obj = {"hostname": mac_based_id}
        entry_point = address

        print(f"Checking registration at: {entry_point}")

        # Make the POST request
        response = requests.post(
            entry_point,
            headers={"Content-Type": "application/json"},
            json=json_obj,
            timeout=100,
        )

        print(f"Registration check response: {response.status_code}")

        if response.status_code == 200:
            try:
                response_data = response.json()
            except json.JSONDecodeError:
                print("Invalid JSON response from server.")
                # myGUIDATAwarn("Invalid JSON response from server.")
                return False

            status = response_data.get("status", False)

            if isinstance(status, str) and status.strip().lower() == "unregistered":
                print("Machine is unregistered. Opening registration page...")
                registration_url = response_data.get("registration_url", "")
                if registration_url:
                    webbrowser.open(registration_url)
                else:
                    print("No registration URL provided.")
                # myGUIDATAwarn("No registration URL provided.")
                return False

            elif isinstance(status, str) and status.strip().lower() == "registered":
                print("Machine is registered. Proceeding...")
                return True

            elif isinstance(status, bool) and not status:
                print("Registration status unclear.")
                # myGUIDATAwarn("Registration status unclear.")
                return False

        else:
            print(
                f"Registration check failed: {response.status_code} - {response.text}"
            )
            # myGUIDATAwarn(f"Registration check failed: {response.status_code} - {response.text}")

    except requests.RequestException as e:
        print(f"Network error during registration check: {e}")
        print("Proceeding without registration check...")
        # myGUIDATAwarn("Network error during registration check. Proceeding without registration check.")
        return True  # Allow offline usage
    except Exception as e:
        print(f"Error during registration check: {e}")
        # myGUIDATAwarn(f"Error during registration check: {e}")

    return False


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


if __name__ == "__main__":

    log_launch_diagnostics()

    # local_remote = "http://127.0.0.1:5000"
    local_remote = "https://test-simplenmr.pythonanywhere.com"
    # local_remote = "https://simplenmr.pythonanywhere.com"

    # At the start of your main program
    app = init_qt_app()

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

    jeol_dict = jeol_data.createJsonDict()

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
            window = MainWindow(str(submission.html_path))
            window.show()
            # Blocks until the viewer window is closed, then ends the program.
            sys.exit(app.exec())
        elif submission.outcome == SubmissionOutcome.DIAGNOSTIC_HTML:
            print(
                f"Server returned a diagnostic report (status {submission.status_code}). "
                "Opening it in the viewer..."
            )
            window = MainWindow(
                str(submission.html_path),
                status_note=f"Server status {submission.status_code} — diagnostic report, not a completed result",
            )
            window.show()
            sys.exit(app.exec())
        else:
            # submit_to_server() has already shown the appropriate dialog
            # (network error / server error / registration required or
            # expired) - nothing further to do here except exit. The JSON
            # file was already saved locally regardless.
            print(f"Server submission did not succeed: {submission.outcome.value}")
            sys.exit(1)
