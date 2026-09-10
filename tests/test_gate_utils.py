import json

from gate_utils import GateLogger, write_gate_status, enforce_strict


def test_gate_logger_prints_and_records(capsys):
    gl = GateLogger()
    gl.log("hello")
    gl.log("world")
    captured = capsys.readouterr()
    assert "hello\nworld" in captured.out
    assert gl.lines == ["hello", "world"]


def test_gate_logger_write_report(tmp_path):
    gl = GateLogger()
    gl.log("line one")
    gl.log("line two")
    out = tmp_path / "report.txt"
    gl.write_report(out)
    assert out.read_text() == "line one\nline two"


def test_write_gate_status_round_trips(tmp_path):
    out = tmp_path / "gate_status.json"
    write_gate_status(out, "PASS", {"system": "wt", "fmax_kj_per_mol_nm": 823.4})
    payload = json.loads(out.read_text())
    assert payload["gate"] == "PASS"
    assert payload["system"] == "wt"
    assert payload["fmax_kj_per_mol_nm"] == 823.4


def test_enforce_strict_exits_on_fail():
    import pytest
    with pytest.raises(SystemExit):
        enforce_strict(True, "FAIL", "boom")


def test_enforce_strict_does_not_exit_on_pass():
    enforce_strict(True, "PASS", "unreachable")


def test_enforce_strict_does_not_exit_when_not_strict():
    enforce_strict(False, "FAIL", "unreachable")
