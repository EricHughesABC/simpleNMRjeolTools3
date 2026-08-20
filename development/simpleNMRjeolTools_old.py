from pathlib import Path
import sys
import uuid
from typing import Dict, List
import requests
import webbrowser
import threading
import h5py
import numpy as np
import pandas as pd
import json
from rdkit import Chem
import fire

from guidata.qthelpers import qt_app_context

from chooseSpectra import NMRExperimentDialog

from qtpy.QtWidgets import QProgressDialog, QApplication, QMessageBox, QFileDialog, QDialog
from qtpy.QtCore import Qt



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


def commandline(fn=None):
    """
    Handle command line arguments for JEOL NMR file processing

    Args:
        fn: Path to a JEOL .jjh5 NMR file (optional)
    """
    # Ensure we have a Qt application context
    with qt_app_context():
        # Convert to Path object if provided
        if fn:
            fn = Path(fn)

        # Validate the provided file
        is_valid, message = validate_file(fn)

        if not is_valid:
            print("Invalid input: opening file dialog...")

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
    """Create the c13predictions dictionary from the c13shifts DataFrame."""
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

    mac_based_id = hex(uuid.getnode())

    hostname = {"datatype": "hostname", "count": 1, "data": {"0": mac_based_id}}

    return hostname


# function to find datasets with peaks
def find_datasets_with_peaks(file_path: Path) -> list[str]:
    """Find datasets with peaks in the HDF5 file.

    Args:
        file_path (Path): Path to the HDF5 file.

    Returns:
        list[str]: List of dataset IDs with peaks.
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
    """Find experiment details for the given datasets.

    Args:
        fn (Path): Path to the HDF5 file.
        expt_id (str): Experiment ID.

    Returns:
        dict: Dictionary of experiment details.
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

            # datasets[expt_id]["specInfo"] = specInfo
        except KeyError:
            print(f"KeyError: SpecInfo not found for experiment {expt_id}")

    return specInfo


def get_peaklist(filename: Path, expt_id: str) -> list:
    try:
        with h5py.File(filename, "r") as fp:
            return list(fp[f"JasonDocument/NMR/NMRData/{expt_id}/Peaks/PeakList"])
    except KeyError:
        return []


def get_peakinfo(filename: Path, expt_id: str, pk_id: str) -> dict:
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


def createWorkingDirectory(fn):
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


def createWorkingFilename(fn):
    # Create the working filename structure
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


# create carbonAtomsinfo output for json
def create_carbonatomsinfo(atoms_df):
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


# create allatomsinfo output for json
def create_allatomsinfo(atoms_df):
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


def readJEOLmolecule(fn):
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
                atom_list.append(atom_info)
                # add atom_idx to atom_list
                atom_info["atom_idx"] = int(atom_idx)
            atoms_df = pd.DataFrame(atom_list)
            # add Z axis all 0
            atoms_df["Z"] = 0.0
        except KeyError as e:
            print(f"KeyError: {e}")

    atoms_df = atoms_df.sort_values(by="atom_idx").reset_index(drop=True)
    return atoms_df





def createMNOVAcalcMethod():
    # Create the MNOVA calculation method structure
    # "MNOVA Predict"
    # "NMRSHIFTDB2 Predict"
    # JEOL Predict
    mNOVA_calc_method = {
        "datatype": "MNOVAcalcMethod",
        "count": 1,
        "data": {"0": "JEOL Predict"},
    }
    return mNOVA_calc_method


def createCarbonCalcPositionsMethod():
    # Create the carbon calculation positions method structure
    carbon_calc_positions_method = {
        "datatype": "carbonCalcPositionsMethod",
        "count": 1,
        "data": {"0": "Calculated Positions"},
    }
    return carbon_calc_positions_method


def dataframe_to_rdkit_molecule(df):
    """
    Convert a pandas DataFrame with molecular structure data directly to RDKit molecule

    Args:
        df: DataFrame with columns ['El', 'NB.Conn', 'NB.Num', 'X', 'Y', 'nH', 'Z']

    Returns:
        rdkit.Chem.Mol: RDKit molecule object or None if creation failed
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
        atom_indices = {}
        for idx, row in mol_df.iterrows():
            element_num = int(row["El"])

            # Create atom with the correct element
            atom = Chem.Atom(element_num)

            # Add atom to molecule
            atom_idx = mol.AddAtom(atom)
            atom_indices[idx] = atom_idx

            element_symbol = atom.GetSymbol()

        # Step 2: Add bonds based on connectivity
        bonds_added = set()  # Track bonds to avoid duplicates
        bond_count = 0

        for idx, row in mol_df.iterrows():
            rdkit_atom_idx = atom_indices[idx]

            # Get connection data
            nb_num = row.get("NB.Num")
            nb_conn = row.get("NB.Conn")

            # Check if connection data exists
            nb_num_valid = nb_num is not None and not (
                isinstance(nb_num, float) and pd.isna(nb_num)
            )
            nb_conn_valid = nb_conn is not None and not (
                isinstance(nb_conn, float) and pd.isna(nb_conn)
            )

            if nb_num_valid and nb_conn_valid:
                try:
                    # Handle numpy arrays and ensure we have lists
                    if isinstance(nb_num, np.ndarray):
                        nb_num = nb_num.tolist()
                    elif not isinstance(nb_num, list):
                        nb_num = [nb_num]

                    if isinstance(nb_conn, np.ndarray):
                        nb_conn = nb_conn.tolist()
                    elif not isinstance(nb_conn, list):
                        nb_conn = [nb_conn]

                    # Process each connection
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

                                # Skip self-connections
                                if connected_df_idx == idx:
                                    continue

                                # Check if connected atom exists in our dataframe
                                if connected_df_idx not in atom_indices:
                                    print(
                                        f"   Warning: Connected atom {connected_df_idx} not found in dataframe"
                                    )
                                    continue

                                connected_rdkit_idx = atom_indices[connected_df_idx]

                                # Avoid duplicate bonds
                                bond_pair = tuple(
                                    sorted([rdkit_atom_idx, connected_rdkit_idx])
                                )
                                if bond_pair not in bonds_added:

                                    # Get bond order
                                    if i < len(nb_conn):
                                        bond_order = nb_conn[i]
                                    else:
                                        bond_order = 1

                                    # Convert bond order to RDKit bond type
                                    if bond_order == 1:
                                        rdkit_bond_type = Chem.BondType.SINGLE
                                        bond_name = "SINGLE"
                                    elif bond_order == 2:
                                        rdkit_bond_type = Chem.BondType.DOUBLE
                                        bond_name = "DOUBLE"
                                    elif bond_order == 3:
                                        rdkit_bond_type = Chem.BondType.TRIPLE
                                        bond_name = "TRIPLE"
                                    elif bond_order == 513:  # Aromatic
                                        rdkit_bond_type = Chem.BondType.AROMATIC
                                        bond_name = "AROMATIC"
                                    elif bond_order == 514:  # Double bond
                                        rdkit_bond_type = Chem.BondType.DOUBLE
                                        bond_name = "DOUBLE (from 514)"
                                    else:
                                        # Default to single bond for unknown types
                                        rdkit_bond_type = Chem.BondType.SINGLE
                                        bond_name = (
                                            f"SINGLE (from unknown {bond_order})"
                                        )

                                    # Add the bond
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

        # Add 2D coordinates (ignoring Z since it's typically 0 for 2D molecules)
        if all(col in mol_df.columns for col in ["X", "Y"]):
            conf = Chem.Conformer(len(mol_df))
            conf.Set3D(False)  # Mark this as a 2D conformer

            for idx, row in mol_df.iterrows():
                rdkit_idx = atom_indices[idx]
                x, y = float(row["X"]), float(row["Y"])
                # For 2D molecules, set Z to 0
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

            # Try sanitization with specific operations
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


def create_molfile(molstr):

    molfile = {"datatype": "molfile", "count": 1, "data": {"0": molstr}}
    return molfile


def create_smiles(smiles_str):
    smiles = {"datatype": "smiles", "count": 1, "data": {"0": smiles_str}}
    return smiles


def get_chosenSpectra(specAssignments):

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


def get_spectraWithPeaks(specAssignments):

    spectraWithPeaks = {
        "datatype": "spectraWithPeaks",
        "count": len(specAssignments),
        "data": {},
    }
    for i, data in enumerate(specAssignments):
        spectraWithPeaks["data"][str(i)] = data["experiment_name"]

    return spectraWithPeaks


class jeolData:
    def __init__(self, file_path: Path):

        self.analysis_cancelled = False
        self.spectra_assignments = {}
        self.file_path = file_path

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
            data["multiplets"]["datatype"] = "multiplets"

        # add workingDirectory
        self.workingDirectory = createWorkingDirectory(self.file_path)

        # add workingFilename
        self.workingFilename = createWorkingFilename(self.file_path)

        # add MNOVAcalcMethod
        self.MNOVAcalcMethod = createMNOVAcalcMethod()

        # add carbonCalcPositionsMethod
        self.carbonCalcPositionsMethod = createCarbonCalcPositionsMethod()

    def updateDatawithExptNames(self):
        for expt_id, data in self.datasets_with_peaks.items():
            print(f"Experiment ID: {data['spec_info']['experimenttype']}")
            # datafilename
            print(f"Data Filename: {data['spec_info']['datafilename']}")
            # 'nucleus'
            print(f"Nucleus: {data['spec_info']['nucleus']}")
            # pulseSequence
            print(f"Pulse Sequence: {data['spec_info']['pulsesequence']}")
            print()

    def createJsonDict(self):
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

    def append_specdata(self):

        expts_type_count = {
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

        expt_kys = []
        expt_types = []
        for idx, expt_id in jeol_data.chosenSpectra["data"].items():
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

        exptNameToKey = {}

        for id, data in jeol_data.datasets_with_peaks.items():
            exptNameToKey[data["spec_info"]["expt_fn"]] = id

        spectra = {}

        for expt_fn, expt_type in zip(expt_kys, expt_types):
            data = None
            for id, dataset in jeol_data.datasets_with_peaks.items():
                if dataset["spec_info"]["expt_fn"] == expt_fn:
                    expt_ky = id
                    data = dataset
                    break
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

        # Use guidata's Qt application context
        with qt_app_context():

            # Generate random experiment names (simulating found experiments)
            experiment_names = []
            chosen_types = {}

            # create experiment_names infor from self.datasets_with_peaks information

            for id, data in self.datasets_with_peaks.items():

                expt_type = data["spec_info"]["experimenttype"]
                expt_fn = Path(data["spec_info"]["datafilename"]).name
                expt_nucleus = data["spec_info"]["nucleus"].__str__()
                # remove quotes from expt_nucleus
                expt_nucleus = expt_nucleus.replace('"', "")
                expt_nucleus = expt_nucleus.replace("'", "")

                expt_name = (
                    f'{expt_nucleus} {data["spec_info"]["pulsesequence"]} {expt_fn}'
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
                    print(
                        f"'{assignment['experiment_name']}' → {assignment['experiment_type']}"
                    )

                # create chosenSpectra data
                self.chosenSpectra = get_chosenSpectra(self.spectra_assignments)

                # create specrtraWithPeaks data
                self.spectraWithPeaks = get_spectraWithPeaks(self.spectra_assignments)

                self.simulated_annealing, self.ml_consent = (
                    dialog.get_processing_options()
                )

            else:
                print("\nDialog was cancelled by user.")
                self.analysis_cancelled = True

            print("\nProgram ending.")


def submit_to_server(json_data: Dict, simpleNMR_address: str) -> bool:
    """
    Submit the JSON data to the processing server with progress dialog.

    Args:
        json_data: The converted JSON data

    Returns:
        True if successful, False otherwise
    """

    # Create progress dialog
    progress = QProgressDialog(
        "Submitting data to simpleNMR server...", "Cancel", 0, 400
    )
    progress.setWindowModality(Qt.WindowModal)
    progress.setMinimumDuration(0)
    progress.setCancelButton(None)  # Remove cancel button since we can't easily cancel the request
    progress.show()

    # Variables to store result
    result = {"success": False, "error": None, "finished": False}

    def make_request():
        try:
            print("Submitting data to simpleNMR server...")

            response = requests.post(
                simpleNMR_address,
                headers={"Content-Type": "application/json"},
                json=json_data,
                timeout=100,
            )

            if response.status_code == 200:
                result["success"] = True
            else:
                error_msg = f"Server error: {response.status_code} - {response.text}"
                print(error_msg)
                result["error"] = error_msg

            # replace dummy_title in response.txt with working_filename from json_data
            workingFilename = json_data["workingFilename"]["data"].get(
                "0", "nmr_analysis_result"
            )
            response_text = response.text
            response_text = response_text.replace("dummy_title", workingFilename)

            # Save response to file
            fn_str = json_data["workingDirectory"]["data"].get("0", ".")
            fn_path = Path(fn_str, "html")

            if not fn_path.exists():
                fn_path.mkdir(parents=True, exist_ok=True)

            # add filename to path
            fn_path = Path(fn_path, workingFilename + ".html")

            with open(fn_path, "w", encoding="utf-8") as f:
                f.write(response_text)

            # Open in browser
            webbrowser.open(f"file://{fn_path}")

        except requests.RequestException as e:
            error_msg = f"Network error: {e}"
            print(error_msg)
            result["error"] = error_msg
        except Exception as e:
            error_msg = f"Error submitting to simpleNMR server: {e}"
            print(error_msg)
            result["error"] = error_msg
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

    # Handle the result
    if result["error"]:
        # Show error dialog
        QMessageBox.critical(
            None,
            "Submission Error",
            f"Failed to submit data to server:\n{result['error']}",
        )
        return False

    return result["success"]


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


if __name__ == "__main__":

    local_remote = "http://127.0.0.1:5000"
    # local_remote = "https://test-simplenmr.pythonanywhere.com"

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
    with open("jeol_data.json", "w") as json_file:
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
        if submit_to_server(jeol_dict, simpleNMR_address):

            print("Analysis complete! Check the opened browser window for results.")
        else:
            print("Server submission failed, but JSON file was saved locally.")
            show_info_message("Submission Failed", "Server submission failed")
