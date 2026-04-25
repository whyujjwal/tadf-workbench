# tests/test_dihedral_3d_panel.py
import pytest
pytest.importorskip("PyQt5")
pytest.importorskip("pyqtgraph.opengl")
from tadf_workbench.ui.dihedral_3d_panel import Dihedral3DPanel
from tadf_workbench.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {
        10: (3.4, 2.9), 30: (3.0, 2.95), 60: (3.1, 2.95), 90: (3.4, 3.0),
    })))


def test_slider_snaps_to_nearest_scanned_angle(qtbot, scan):
    panel = Dihedral3DPanel()
    qtbot.addWidget(panel)
    panel.set_scan(scan)
    panel.set_angle(35)  # closest = 30
    assert panel.current_angle == 30.0
    panel.set_angle(75)  # closest = 60 or 90 → 60 vs 90 — 75-60=15 < 90-75=15 tie → parser chooses 60 by min()
    assert panel.current_angle in (60.0, 90.0)
    panel.set_angle(89)
    assert panel.current_angle == 90.0


def test_panel_emits_angle_changed(qtbot, scan):
    panel = Dihedral3DPanel()
    qtbot.addWidget(panel)
    panel.set_scan(scan)
    received = []
    panel.angle_changed.connect(received.append)
    panel.set_angle(31)  # snaps to 30
    assert received and received[-1] == 30.0


def test_setting_same_angle_twice_does_not_double_emit(qtbot, scan):
    panel = Dihedral3DPanel()
    qtbot.addWidget(panel)
    panel.set_scan(scan)
    received = []
    panel.angle_changed.connect(received.append)
    panel.set_angle(30)
    panel.set_angle(30)
    assert received.count(30.0) == 1
