# tests/test_angle_scan_window.py
import pytest
pytest.importorskip("PyQt5")
from tadf_workbench.ui import AngleScanWindow
from tadf_workbench.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {
        10: (3.4, 2.9), 30: (3.0, 2.95), 60: (3.1, 2.95), 90: (3.4, 3.0),
    })))


def test_window_loads_scan(qtbot, scan):
    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)
    assert w.windowTitle().startswith("Angle Scan")
    assert w.table.rowCount() == 4


def test_selecting_angle_in_table_updates_3d(qtbot, scan):
    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)
    # Programmatically select row for θ=30
    angles = [p.angle_deg for p in scan.sorted_points]
    row = angles.index(30.0)
    w.table.selectRow(row)
    qtbot.wait(50)
    assert w.viewer_3d.current_angle == 30.0


def test_selecting_angle_in_3d_updates_plots(qtbot, scan):
    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)
    # Slider change → plots should receive new highlight
    w.viewer_3d.set_angle(60)
    qtbot.wait(50)
    assert w.viewer_3d.current_angle == 60.0
    assert w._current_angle == 60.0


def test_status_bar_reports_min_gap(qtbot, scan):
    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)
    txt = w.statusBar().currentMessage()
    assert "min gap" in txt.lower()
    assert "30" in txt  # min gap is at 30°
