"""
JEOL universal IDs carried into the simpleNMR input JSON.

Builds a tiny HDF5 file shaped like a .jjh5 (one spectrum with SpecInfo,
one peak and one integral, each with an "ID" attribute) and checks that
the readers in json_converter return those IDs.
"""

import h5py
import numpy as np
import pytest

from simplenmrjeol.json_converter import get_integralinfo, get_peakinfo, get_spec_info

SPECTRUM_ID = "{23031bd1-9838-48c4-9280-9a37fc41e3cb}"
PEAK_ID = "{1774e688-a307-4f0c-92ad-b60b6da1c87b}"
INTEGRAL_ID = "{bb4861a6-78c2-48fa-930f-913298c250d2}"


@pytest.fixture
def jjh5(tmp_path):
    path = tmp_path / "mini.jjh5"
    with h5py.File(path, "w") as f:
        spec = f.create_group("JasonDocument/NMR/NMRData/4")
        spec.attrs["ID"] = SPECTRUM_ID.encode()          # JASON stores IDs as bytes

        info = spec.create_group("SpecInfo")
        info.attrs["PulseProgram"] = b"hmbcetgpl3nd"
        info.attrs["ExperimentType.str"] = b"HMBC"
        info.attrs["OrigFilename.filepath.str"] = b"/data/ethylbenzene/16/pdata/1"
        info.attrs["Solvent"] = b"CDCl3"
        info.attrs["SpectrometerFrequencies"] = np.array([598.91, 150.596, 0.0])
        for n, (name, isotope) in enumerate([(b"H", 1), (b"C", 13)]):
            nuc = info.create_group(f"Nucleides/{n}")
            nuc.attrs["Name"] = np.bytes_(name)   # fixed-length bytes, as JASON writes it
            nuc.attrs["Isotope"] = isotope

        peak = spec.create_group("Peaks/PeakList/0")
        peak.attrs["ID"] = PEAK_ID.encode()
        peak.attrs["Height"] = 1.0e6
        peak.attrs["Pos"] = np.array([7.4039, 144.4006, 0.0])

        integral = spec.create_group("Multiplets_Integrals/MultipletList/0")
        integral.attrs["ID"] = INTEGRAL_ID.encode()
        integral.attrs["Value"] = 1.0
        integral.attrs["SpectrumRange[0]"] = np.array([7.31, 7.33])
        integral.attrs["SpectrumRange[1]"] = np.array([127.78, 128.10])
    return path


def test_peak_has_jason_id(jjh5):
    peak = get_peakinfo(jjh5, "4", "0")
    assert peak["jason_id"] == PEAK_ID
    assert (peak["delta2"], peak["delta1"]) == pytest.approx((7.4039, 144.4006))


def test_integral_has_jason_id(jjh5):
    assert get_integralinfo(jjh5, "4", "0")["jason_id"] == INTEGRAL_ID


def test_spectrum_has_jason_ids(jjh5):
    info = get_spec_info(jjh5, "4")
    assert info["jason_spectrum_id"] == SPECTRUM_ID
    assert info["jason_nmrdata_index"] == "4"
    assert info["experimenttype"] == "HMBC"


def test_missing_id_gives_none(jjh5):
    with h5py.File(jjh5, "a") as f:
        del f["JasonDocument/NMR/NMRData/4/Peaks/PeakList/0"].attrs["ID"]
    assert get_peakinfo(jjh5, "4", "0")["jason_id"] is None
