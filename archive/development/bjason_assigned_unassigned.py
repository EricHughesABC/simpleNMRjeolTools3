from pathlib import Path
import platform
import h5py
import beautifuljason as bjason

import numpy as np
import pandas as pd

if __name__ == "__main__":

    fn_assigned = Path("/Users/vsmw51/Downloads/alpha_ionone_assigned_by_jason.jjh5")
    fn_unassigned = Path("/Users/vsmw51/Downloads/alpha_ionone_unassigned.jjh5")

    print(f"fn_assigned exists: {fn_assigned.exists()}")
    print(f"fn_unassigned exists: {fn_unassigned.exists()}")

    with bjason.Document(fn_assigned, mode="r") as doc:
        mol = doc.mol_data[0]
        spec = mol.spectra[0]

        for shift in spec.shifts:
            value, _ = shift.get_value_error_pair(
                bjason.Molecule.CalcMethod.Experimental
            )
            if value is not None:
                print(list(shift.nums), shift.mark, value)