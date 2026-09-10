"""
05_analysis.py's MDAnalysis-dependent run_analysis() and the fpocket
subprocess wrapper aren't covered here -- they're isolated behind local
imports / shutil.which() checks specifically so the pure gate/geometry/
parsing functions stay testable without MDAnalysis or fpocket installed
(neither is in this test environment). These tests cover the parts with
actual logic: trajectory_sanity(), na_coordination_distance(),
parse_fpocket_info(), and select_largest_pocket().
"""
import numpy as np

from conftest import load_script

SAMPLE_FPOCKET_INFO = """\
Pocket 1 :
	Score : 	0.523
	Druggability Score : 	0.812
	Number of Alpha Spheres : 	42
	Total SASA : 	310.221
	Volume : 	850.331
	Charge Score : 	2
	Local hydrophobic density Score : 	28.5

Pocket 2 :
	Score : 	0.201
	Druggability Score : 	0.334
	Number of Alpha Spheres : 	15
	Total SASA : 	120.004
	Volume : 	210.117
	Charge Score : 	0
	Local hydrophobic density Score : 	10.2
"""


def test_trajectory_sanity_zero_frames():
    mod = load_script("05_analysis.py")
    ok, reason = mod.trajectory_sanity(0)
    assert ok is False
    assert "zero frames" in reason


def test_trajectory_sanity_below_minimum():
    mod = load_script("05_analysis.py")
    ok, reason = mod.trajectory_sanity(3, min_frames=10)
    assert ok is False
    assert "3" in reason


def test_trajectory_sanity_enough_frames():
    mod = load_script("05_analysis.py")
    ok, reason = mod.trajectory_sanity(500, min_frames=10)
    assert ok is True


def test_na_coordination_distance():
    mod = load_script("05_analysis.py")
    na_pos = np.array([3.0, 0.0, 0.0])
    target_pos = np.array([0.0, 0.0, 0.0])
    assert mod.na_coordination_distance(na_pos, target_pos) == 3.0


def test_resid438_atom_default_differs_ser_vs_thr():
    """The whole reason this default is a dict keyed by system rather than
    one shared constant: Ser438 (WT) uses OG, Thr438 (S438T) uses OG1.
    Sharing one default would silently find nothing for one system."""
    mod = load_script("05_analysis.py")
    assert mod.RESID_438_ATOM_DEFAULT["wt"] == "OG"
    assert mod.RESID_438_ATOM_DEFAULT["s438t"] == "OG1"
    assert mod.RESID_438_ATOM_DEFAULT["wt"] != mod.RESID_438_ATOM_DEFAULT["s438t"]


def test_parse_fpocket_info_two_pockets():
    mod = load_script("05_analysis.py")
    pockets = mod.parse_fpocket_info(SAMPLE_FPOCKET_INFO)

    assert len(pockets) == 2
    assert pockets[0]["pocket_id"] == 1
    assert pockets[0]["volume_A3"] == 850.331
    assert pockets[0]["score"] == 0.523
    assert pockets[1]["pocket_id"] == 2
    assert pockets[1]["volume_A3"] == 210.117


def test_parse_fpocket_info_empty_text():
    mod = load_script("05_analysis.py")
    assert mod.parse_fpocket_info("") == []


def test_select_largest_pocket_picks_max_volume():
    mod = load_script("05_analysis.py")
    pockets = mod.parse_fpocket_info(SAMPLE_FPOCKET_INFO)

    largest = mod.select_largest_pocket(pockets)

    assert largest["pocket_id"] == 1
    assert largest["volume_A3"] == 850.331


def test_select_largest_pocket_empty_list_returns_none():
    mod = load_script("05_analysis.py")
    assert mod.select_largest_pocket([]) is None
