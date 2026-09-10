"""
Covers the actual bug class this pipeline is trying to prevent: a
CHARMM-GUI/CGenFF export that's missing something (a topology, an index,
a coordinate file, all .mdp files, or the CGenFF ligand parameters
escitalopram needs) should FAIL loudly here, not pass through and break
02_minimize.py three steps later with a confusing gmx error about an
unresolved residue.
"""
from conftest import load_script


def _make_complete_export(root):
    (root / "topol.top").write_text("; dummy topology\n")
    (root / "index.ndx").write_text("[ System ]\n1 2 3\n")
    (root / "step5_input.gro").write_text("dummy\n1\n1.0 1.0 1.0\n")
    (root / "step6.1_equilibration.mdp").write_text("integrator = md\n")
    (root / "step7_production.mdp").write_text("integrator = md\n")
    (root / "escitalopram.str").write_text("* CGenFF stream file\n")


def test_validate_complete_export_passes(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    _make_complete_export(tmp_path)

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["missing"] == []
    assert found["topol"] is not None
    assert found["index"] is not None
    assert found["coordinates"] is not None
    assert len(found["mdp_files"]) == 2
    assert len(found["cgenff_files"]) == 1


def test_validate_missing_topology_fails(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    _make_complete_export(tmp_path)
    (tmp_path / "topol.top").unlink()

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["topol"] is None
    assert any("topology" in m for m in found["missing"])


def test_validate_missing_mdp_files_fails(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    _make_complete_export(tmp_path)
    for mdp in tmp_path.glob("*.mdp"):
        mdp.unlink()

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["mdp_files"] == []
    assert any("mdp" in m for m in found["missing"])


def test_validate_missing_cgenff_params_fails(tmp_path):
    """This is the SERT-specific check psilocin-5ht2a's version of this
    script doesn't need -- escitalopram is a CGenFF ligand, and a download
    missing the .str stream file (or split .rtf/.prm) will resolve
    everything else fine and then fail deep inside grompp on an unknown
    residue, which is a much worse place to find out."""
    mod = load_script("01_import_charmm_gui.py")
    _make_complete_export(tmp_path)
    (tmp_path / "escitalopram.str").unlink()

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["cgenff_files"] == []
    assert any("CGenFF" in m for m in found["missing"])


def test_validate_accepts_split_rtf_prm_instead_of_stream(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    _make_complete_export(tmp_path)
    (tmp_path / "escitalopram.str").unlink()
    (tmp_path / "escitalopram.rtf").write_text("* dummy rtf\n")
    (tmp_path / "escitalopram.prm").write_text("* dummy prm\n")

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["missing"] == []
    assert len(found["cgenff_files"]) == 2


def test_validate_accepts_alternate_coordinate_name(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    (tmp_path / "topol.top").write_text("; dummy\n")
    (tmp_path / "index.ndx").write_text("[ System ]\n")
    (tmp_path / "system.gro").write_text("dummy\n")
    (tmp_path / "step7_production.mdp").write_text("integrator = md\n")
    (tmp_path / "escitalopram.str").write_text("* dummy\n")

    found = mod.validate_charmm_gui_dir(tmp_path)

    assert found["coordinates"] is not None
    assert found["missing"] == []


def test_stage_into_copies_found_files(tmp_path):
    mod = load_script("01_import_charmm_gui.py")
    src = tmp_path / "src"
    src.mkdir()
    _make_complete_export(src)
    found = mod.validate_charmm_gui_dir(src)

    out_dir = tmp_path / "staged"
    mod.stage_into(out_dir, found)

    assert (out_dir / "topol.top").exists()
    assert (out_dir / "index.ndx").exists()
    assert (out_dir / "step5_input.gro").exists()
    assert (out_dir / "step6.1_equilibration.mdp").exists()
    assert (out_dir / "step7_production.mdp").exists()
    assert (out_dir / "escitalopram.str").exists()
