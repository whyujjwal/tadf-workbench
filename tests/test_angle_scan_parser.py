# tests/test_angle_scan_parser.py
import pytest
from tadf_workbench.parsers.angle_scan_parser import angle_from_filename


@pytest.mark.parametrize("name,expected", [
    ("10.log", 10.0),
    ("angle_30.log", 30.0),
    ("scan_d45deg.out", 45.0),
    ("mol_60_opt.log", 60.0),
    ("90deg.log", 90.0),
    ("theta-75.log", 75.0),
    ("phi_0_opt.out", 0.0),
    ("donor_acceptor_scan_180.log", 180.0),
])
def test_angle_from_filename_extracts_first_integer(name, expected):
    assert angle_from_filename(name) == expected


def test_angle_from_filename_handles_no_number():
    assert angle_from_filename("no_number_here.log") is None
    assert angle_from_filename("final.log") is None


def test_angle_from_filename_ignores_directory_numbers():
    assert angle_from_filename("/data/scan_2025/30.log") == 30.0
    assert angle_from_filename("/2024/final/angle_45.log") == 45.0


from pathlib import Path
from tadf_workbench.parsers import parse_angle_folder
from tadf_workbench.core import AngleScan
from tests.fixtures.build_fixtures import write_scan_folder


def test_parse_angle_folder_reads_all_logs(tmp_path):
    folder = write_scan_folder(tmp_path / "scan", {
        10: (3.40, 2.90),
        30: (3.35, 2.90),
        60: (3.20, 3.10),
    })
    scan = parse_angle_folder(str(folder))
    assert isinstance(scan, AngleScan)
    assert len(scan) == 3
    assert scan.angles == [10.0, 30.0, 60.0]
    p30 = next(p for p in scan.points if p.angle_deg == 30)
    assert p30.s1_t1_gap_ev == pytest.approx(0.45, abs=1e-3)


def test_parse_angle_folder_skips_non_log_files(tmp_path):
    folder = tmp_path / "mixed"
    write_scan_folder(folder, {10: (3.4, 2.9)})
    (folder / "notes.txt").write_text("ignore me")
    (folder / "README").write_text("ignore me too")
    scan = parse_angle_folder(str(folder))
    assert len(scan) == 1


def test_parse_angle_folder_skips_files_without_angle(tmp_path):
    folder = tmp_path / "scan"
    folder.mkdir()
    # Will contain a valid log body but filename has no number
    from tests.fixtures.build_fixtures import write_fixture
    write_fixture(folder / "final_opt.log", s1_energy_ev=3.4, t1_energy_ev=2.9)
    write_fixture(folder / "angle_30.log", s1_energy_ev=3.3, t1_energy_ev=2.9)
    scan = parse_angle_folder(str(folder))
    assert len(scan) == 1
    assert scan.points[0].angle_deg == 30.0


def test_parse_angle_folder_raises_on_missing_dir(tmp_path):
    with pytest.raises(FileNotFoundError):
        parse_angle_folder(str(tmp_path / "does_not_exist"))
