from abc import ABC, abstractmethod
import pandas as pd
import uuid
from pathlib import Path
import h5py
from rdkit import Chem

class BaseNMRExtractor(ABC):
    def __init__(self, file_path):
        self.file_path = file_path


    @abstractmethod
    def setup(self, fn: Path) -> None:
        """
        Perform any necessary setup or initialization before data extraction.
        This method can be used to prepare the environment, load necessary libraries,
        or perform any pre-processing steps required for the specific NMR software.

        Args:
            fn (Path): The file path to the NMR data file that will be processed.
        """
        pass

    # --- 1. File & Environment Metadata ---
    @abstractmethod
    def hostname(self) -> dict: pass

    @abstractmethod
    def workingDirectory(self, fn: Path) -> dict: pass

    @abstractmethod
    def workingFilename(self, fn: Path) -> dict: pass

    # --- 2. Chemical Structure & Atoms ---
    @abstractmethod
    def smiles(self, fn: Path) -> dict: pass

    @abstractmethod
    def molfile(self, fn: Path) -> dict: pass

    @abstractmethod
    def allAtomsInfo(self) -> dict: pass

    @abstractmethod
    def carbonAtomsInfo(self) -> dict: pass

    # --- 3. NMR Data & Assignments ---
    @abstractmethod
    def nmrAssignments(self) -> dict: pass

    @abstractmethod
    def c13predictions(self) -> dict: pass

    # --- 4. Spectral Data (1D/2D) ---
    @abstractmethod
    def spectrum(self) -> dict: pass

    @abstractmethod
    def chosenSpectra(self) -> dict: pass

    @abstractmethod
    def exptIdentifiers(self) -> dict: pass

    @abstractmethod
    def spectraWithPeaks(self) -> dict: pass

    # --- 5. Calculation & Simulation Settings ---
    @abstractmethod
    def carbonCalcPositionsMethod(self) -> dict: pass

    @abstractmethod
    def MNOVAcalcMethod(self) -> dict: pass

    @abstractmethod
    def randomizeStart(self) -> dict: pass

    @abstractmethod
    def startingTemperature(self) -> dict: pass

    @abstractmethod
    def endingTemperature(self) -> dict: pass

    @abstractmethod
    def coolingRate(self) -> dict: pass

    @abstractmethod
    def numberOfSteps(self) -> dict: pass

    @abstractmethod
    def ppmGroupSeparation(self) -> dict: pass

    @abstractmethod
    def ml_consent(self) -> dict: pass

    @abstractmethod
    def simulatedAnnealing(self) -> dict: pass

    @abstractmethod
    def append_specdata(self) -> dict: pass

    # # --- Standardized Output Method ---
    # def generate_payload(self) -> dict:
    #     """
    #     This is a concrete method (not abstract). 
    #     It builds the final dictionary to send to the server.
    #     """
    #     return {
    #         "metadata": {
    #             "host": self.hostname(),
    #             "dir": self.workingDirectory(),
    #             "file": self.workingFilename()
    #         },
    #         "structure": {
    #             "smiles": self.smiles(),
    #             "mol": self.molfile(),
    #             "atoms": self.allAtomsInfo(),
    #             "carbons": self.carbonAtomsInfo()
    #         },
    #         "nmr": {
    #             "assignments": self.nmrAssignments(),
    #             "predictions": self.c13predictions(),
    #             "spectra": self.spectrum()
    #         },
    #         "settings": {
    #             "temp_range": (self.startingTemperature(), self.endingTemperature()),
    #             "steps": self.numberOfSteps(),
    #             "ml_consent": self.ml_consent()
    #         }
    #     }
    
    def createSimpleNMRjsonDict(self):
        """
        Creates and returns a dictionary containing all relevant data for JSON serialization.

        The dictionary includes molecular information, calculation methods, atom data, predictions,
        spectra information, simulated annealing results, machine learning consent, and additional
        spectral data appended from `append_specdata()`.

        Returns:
            dict: A dictionary with all necessary fields for JSON output.
        """

        json_dict = {
            "smiles": self.smiles(self.file_path),
            "molfile": self.molfile(self.file_path),
            "hostname": self.hostname(),
            "workingDirectory": self.workingDirectory(self.file_path),
            "workingFilename": self.workingFilename(self.file_path),
            "MNOVAcalcMethod": self.MNOVAcalcMethod(),
            "carbonCalcPositionsMethod": self.carbonCalcPositionsMethod(),
            "allAtomsInfo": self.allAtomsInfo(),
            "carbonAtomsInfo": self.carbonAtomsInfo(),
            "c13predictions": self.c13predictions(),
            "chosenSpectra": self.chosenSpectra(),
            "spectraWithPeaks": self.spectraWithPeaks(),
            "simulatedAnnealing": {
                "datatype": "simulatedAnnealing",
                "count": 1,
                "data": {"0": self.simulatedAnnealing()},
            },
            "ml_consent": {
                "datatype": "ml_consent",
                "count": 1,
                "data": {"0": self.ml_consent()},
            },
        }

        spec_dict = self.append_specdata()

        # merge the two dictionaries
        json_dict.update(spec_dict)

        return json_dict
    
class JeolNMRExtractor(BaseNMRExtractor):
    def __init__(self, file_path):
        super().__init__(file_path)
        # Initialize any Jeol-specific attributes here

    def _create_all_atoms_info_from_mol(self) -> Dict[str, Any]:
        """
        Create the allAtomsInfo structure from the RDKit molecule.
        
        Returns:
            Dictionary containing all atom information
        """
        if not RDKIT_AVAILABLE or not self.rdkit_mol:
            return {"datatype": "allAtomsInfo", "data": {}, "count": 0}
        
        all_atoms_data = {
            "datatype": "allAtomsInfo",
            "data": {},
            "count": self.rdkit_mol.GetNumAtoms()
        }
        
        for atom_idx, atom in enumerate(self.rdkit_mol.GetAtoms()):
            atom_info = {
                "atom_idx": atom_idx,
                "id": atom_idx,
                "atomNumber": str(atom_idx + 1),  # 1-based numbering as string
                "symbol": atom.GetSymbol(),
                "numProtons": atom.GetTotalNumHs()
            }
            all_atoms_data["data"][str(atom_idx)] = atom_info
        
        return all_atoms_data


    def readJEOLmolecule(self, fn: Path) -> pd.DataFrame:
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

    # Implement all abstract methods with Jeol-specific logic
    def hostname(self):
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


    def workingDirectory(self, fn: Path) -> dict:
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

    def workingFilename(self, fn: Path) -> dict:
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

    def smiles(self, fn: Path) -> dict:
        # Extract and return the SMILES string from the Jeol file
        
        molfile_dict = self.molfile(fn)

        if molfile_dict["count"] == 0:
            print("No molfile found, cannot extract SMILES.")
            return {"datatype": "smiles", "count": 0, "data": {}}
        
        molfile_str = molfile_dict["data"]["0"]
        # Here you would implement the actual conversion from molfile to SMILES.
        smiles_str = Chem.MolToSmiles(molfile_dict["data"]["0"])
        return {"datatype": "smiles", "count": 1, "data": {"0": smiles_str}}

    def molfile(self, fn: Path) -> dict:
        # Extract and return the molfile from the Jeol file

        with h5py.File(fn, 'r') as fp:
            try:
                mol_info = fp['JasonDocument/Molecules/Molecules/0']

                molfilestr = mol_info.attrs['sdfile']

            except KeyError:
                print("Key not found")
                return {"datatype": "molfile", "count": 0, "data": {}}
            
        if isinstance(molfilestr, bytes):
            molfilestr = molfilestr.decode('utf-8')

        return {"datatype": "molfile", "count": 1, "data": {"0": molfilestr}}


    def allAtomsInfo(self):
        # Extract and return information about all atoms from the Jeol file
        pass

    def carbonAtomsInfo(self):
        # Extract and return information about carbon atoms from the Jeol file
        pass

    def nmrAssignments(self):
        # Extract and return NMR assignments from the Jeol file
        pass

    def c13predictions(self):
        # Extract and return C13 predictions from the Jeol file
        pass

    def spectrum(self):
        # Extract and return spectrum data from the Jeol file
        pass

    def chosenSpectra(self):
        # Extract and return chosen spectra from the Jeol file
        pass

    def exptIdentifiers(self):
        # Extract and return experiment identifiers from the Jeol file
        pass

    def spectraWithPeaks(self):
        # Extract and return spectra with peaks from the Jeol file
        pass

    def carbonCalcPositionsMethod(self):
        # Return the method used for calculating carbon positions in Jeol
        pass

    def MNOVAcalcMethod(self):
        # Return the method used for MNOVA calculations in Jeol
        pass

    def randomizeStart(self):
        # Return whether randomize start is enabled for Jeol
        pass

    def startingTemperature(self):
        # Return the starting temperature for simulated annealing in Jeol
        pass

    def endingTemperature(self):
        # Return the ending temperature for simulated annealing in Jeol
        pass

    def coolingRate(self):
        # Return the cooling rate for simulated annealing in Jeol
        pass

    def numberOfSteps(self):
        # Return the number of steps for simulated annealing
        pass
