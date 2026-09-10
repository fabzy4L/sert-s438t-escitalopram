import numpy as np

from conftest import load_script


def test_inputs_compatible_matching_resids():
    mod = load_script("06_compare_systems.py")
    ok, reason = mod.inputs_compatible([95, 98, 438, 439], [95, 98, 438, 439])
    assert ok is True


def test_inputs_compatible_mismatched_resids():
    mod = load_script("06_compare_systems.py")
    ok, reason = mod.inputs_compatible([95, 98, 438, 439], [95, 98, 438])
    assert ok is False
    assert "differ" in reason


def test_inputs_compatible_empty_fails():
    mod = load_script("06_compare_systems.py")
    ok, reason = mod.inputs_compatible([], [95, 98, 438, 439])
    assert ok is False


def test_compare_rmsf_matches_by_resid_not_position():
    """The two systems' 05_analysis.py runs aren't guaranteed to list
    residues in the same order -- this must match by resid, not zip by
    array position, or a S438T-vs-WT delta could silently compare the
    wrong pair of residues."""
    mod = load_script("06_compare_systems.py")
    wt_resids, wt_rmsf = [95, 98, 438, 439], np.array([1.0, 1.2, 0.8, 1.1])
    # deliberately different order
    mut_resids, mut_rmsf = [439, 438, 98, 95], np.array([1.3, 1.5, 1.1, 0.9])

    rows = mod.compare_rmsf(wt_resids, wt_rmsf, mut_resids, mut_rmsf)
    by_resid = {r["resid"]: r for r in rows}

    assert by_resid[438]["rmsf_wt_angstrom"] == 0.8
    assert by_resid[438]["rmsf_s438t_angstrom"] == 1.5
    assert abs(by_resid[438]["delta_angstrom"] - 0.7) < 1e-9
    assert by_resid[95]["rmsf_wt_angstrom"] == 1.0
    assert by_resid[95]["rmsf_s438t_angstrom"] == 0.9


def test_compare_na_coordination_computes_delta():
    mod = load_script("06_compare_systems.py")
    wt = np.array([2.5, 2.6, 2.4, 2.5])
    mut = np.array([4.0, 4.2, 3.9, 4.1])

    summary = mod.compare_na_coordination(wt, mut)

    assert abs(summary["wt_mean_angstrom"] - 2.5) < 1e-9
    assert abs(summary["s438t_mean_angstrom"] - 4.05) < 1e-9
    assert summary["delta_mean_angstrom"] > 0
    assert summary["wt_n"] == 4
    assert summary["s438t_n"] == 4


def test_compare_na_coordination_empty_series():
    """Also pins the empty-branch dict to the SAME keys as the populated
    branch (wt_mean_angstrom, not wt_mean) -- these two branches returning
    different key sets is exactly the kind of mismatch that would
    KeyError in main()'s unconditional na_summary["wt_mean_angstrom"]
    lookup, silently, only on the empty-data path."""
    mod = load_script("06_compare_systems.py")
    summary = mod.compare_na_coordination(np.array([]), np.array([1.0]))
    assert summary["wt_mean_angstrom"] is None
    assert summary["s438t_mean_angstrom"] is None
    assert set(summary.keys()) == {
        "wt_mean_angstrom", "wt_std_angstrom", "s438t_mean_angstrom",
        "s438t_std_angstrom", "delta_mean_angstrom", "wt_n", "s438t_n",
    }
