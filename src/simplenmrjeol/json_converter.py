"""
json_converter.py

Reads a JEOL/JASON .jjh5 NMR file (molecule, HSQC/HMBC/COSY/etc. peak
lists, C13 predictions) and builds the JSON payload simpleNMR's Flask
server expects — the same job Bruker's json_converter.py does for
TopSpin data, hence the matching name.

Everything here is data extraction and JSON construction. The actual
runnable program (dialog flow, server submission, opening the results
viewer) lives in jason_simpleNMR_cli.py, which imports jeolData from
this module.
"""

from __future__ import annotations

from pathlib import Path
import uuid

import h5py
import numpy as np
import pandas as pd
from rdkit import Chem

from qtpy.QtWidgets import QDialog

# simplenmr_builder is a hard dependency (see jason_simpleNMR_cli.py for
# the user-facing "not installed" error — this module just imports
# normally and lets ImportError propagate naturally, as any library
# module should; the entry point is responsible for a clear message).
from simplenmr_builder import SimpleNMRBuilder
from simplenmr_builder.gui.spectrum_assignment_dialog import (
    SpectrumAssignmentDialog as NMRExperimentDialog,
)

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

    def append_specdata(self, builder: "SimpleNMRBuilder") -> None:
        """
        Same experiment-matching logic as append_specdata() (splitting
        each chosenSpectra entry, looking up the matching dataset,
        attaching spec_info/peaks/integrals/multiplets), but adds each
        block through builder.spectra.add_block(token, block) instead of
        indexing a hand-maintained expts_type_count dict.

        This is the exact mechanism that caused the confirmed 2026-08-27
        KeyError bug for DEPT135/DDEPTCH3ONLY (see CHANGES.md) —
        SpectrumBuilder owns per-token block numbering internally and
        validates every token against the canonical NMREXPERIMENTS list,
        so an unrecognized token now raises a clear ValueError immediately
        instead of crashing deep in a dict lookup, or (for a token that's
        merely wrong rather than absent from the counter dict) silently
        vanishing with no error anywhere.
        """
        expt_kys: list[str] = []
        expt_types: list[str] = []
        for idx, expt_id in self.chosenSpectra["data"].items():
            expt_id_parts = expt_id.split(" ")
            if len(expt_id_parts) > 4:
                expt_id_parts = [" ".join(expt_id_parts[:2])] + expt_id_parts[2:]
            expt_kys.append(expt_id_parts[-2])
            expt_types.append(expt_id_parts[-1])

        for expt_ky, expt_type in zip(expt_kys, expt_types):
            data = self.datasets_with_peaks.get(expt_ky, None)
            if not data:
                print(f"Experiment {expt_ky} not found")
                continue

            block = {"datatype": "nmrspectrum"}
            for k, v in data["spec_info"].items():
                block[k] = v
            block["peaks"] = data["peaks"]
            block["integrals"] = data["integrals"]
            block["multiplets"] = data["multiplets"]

            builder.spectra.add_block(expt_type, block)

    def createJsonDict(self, validator_path=None) -> dict:
        """
        Preferred replacement for createJsonDict(): builds the submission
        payload through simplenmr_builder.SimpleNMRBuilder instead of
        hand-assembling json_dict directly. Extracts raw values out of the
        envelope dicts already computed in __init__ (self.smiles,
        self.molfile, etc. are already {"datatype","count","data"}
        envelopes from create_smiles()/create_molfile()/etc. — the builder
        wants raw values and rebuilds the envelope itself, so this reads
        ["data"]["0"] or list(...["data"].values()) as appropriate rather
        than re-deriving anything).

        Runs the real required-field checks (refuses to proceed if
        smiles/molfile/carbonAtomsInfo would be missing — the three fields
        with no server-side fallback) and the real schema + HSQC-presence
        validation from validate_simplenmr_json.py before returning,
        matching what the server itself will check.
        """
        builder = SimpleNMRBuilder(source="jeol")

        builder.set_scalar("smiles", self.smiles["data"]["0"])
        builder.set_scalar("molfile", self.molfile["data"]["0"])
        builder.set_scalar("hostname", self.hostname["data"]["0"])
        builder.set_scalar("MNOVAcalcMethod", self.MNOVAcalcMethod["data"]["0"])
        builder.set_scalar("carbonCalcPositionsMethod", self.carbonCalcPositionsMethod["data"]["0"])
        builder.set_scalar("simulatedAnnealing", self.simulated_annealing)
        builder.set_scalar("ml_consent", self.ml_consent)
        # workingDirectory/workingFilename: server_consumed=false per the
        # manifest (harmless passthrough server-side), but submit_to_server()
        # / _write_html_result() / _log_submission_error() below all read
        # these back OUT of json_data locally, after the round trip, to know
        # where to write the HTML result and error log. Omitting them here
        # doesn't break the server call itself but crashes the *local*
        # post-submission file-writing step with a bare KeyError — confirmed
        # 2026-08-27 from a real run ("Error submitting to simpleNMR server:
        # 'workingFilename'"), even though the server had already responded
        # successfully. Both fields must be included regardless of whether
        # the schema marks them required.
        builder.set_scalar("workingDirectory", self.workingDirectory["data"]["0"])
        builder.set_scalar("workingFilename", self.workingFilename["data"]["0"])

        builder.set_all_atoms_info(list(self.allatomsInfo["data"].values()))
        builder.set_carbon_atoms_info(list(self.carbonatomsInfo["data"].values()))
        builder.set_c13predictions(list(self.c13predictions["data"].values()))

        self.append_specdata(builder)

        # self.chosenSpectra was already filtered for SKIP entries by the
        # fix in choosePeakPickedSpectaforSimpleNMR() above (list
        # comprehension, not mutate-while-iterating) — skip=False here
        # just registers each already-filtered entry with the builder for
        # consistency, it doesn't re-filter anything.
        for entry in self.chosenSpectra["data"].values():
            builder.spectra.add_chosen_candidate(entry, skip=False)

        for entry in self.spectraWithPeaks["data"].values():
            builder.spectra.add_spectrum_with_peaks(entry)

        return builder.build(validator_path=validator_path)

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

        dialog = NMRExperimentDialog(experiment_names, chosen_types=chosen_types)

        # Execute dialog and get result
        result = dialog.exec_()

        if result == QDialog.Accepted:
            self.spectra_assignments = dialog.get_experiment_assignments()

            # Filter out SKIP-assigned spectra by building a new list rather
            # than removing from self.spectra_assignments while iterating
            # over it. The previous in-place .remove() approach shifted list
            # indices after every removal, which caused the loop to skip the
            # element immediately following a removed one — so when two or
            # more SKIP-assigned spectra were adjacent in the list, some of
            # them silently survived into chosenSpectra and got submitted as
            # real experiment data despite being explicitly marked SKIP.
            # See simpleNMR_field_manifest.yaml for the full writeup and
            # test_skip_filtering.py for a standalone reproduction.
            self.spectra_assignments = [
                assignment
                for assignment in self.spectra_assignments
                if assignment["experiment_type"] != "SKIP"
            ]

            # create chosenSpectra data
            self.chosenSpectra = get_chosenSpectra(self.spectra_assignments)

            # create specrtraWithPeaks data
            self.spectraWithPeaks = get_spectraWithPeaks(self.spectra_assignments)

            self.simulated_annealing, self.ml_consent = dialog.get_processing_options()

        else:
            print("\nDialog was cancelled by user.")
            self.analysis_cancelled = True
