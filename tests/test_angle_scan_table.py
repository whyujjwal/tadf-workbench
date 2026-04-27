# tests/test_angle_scan_table.py
import pytest
pytest.importorskip("PyQt5")
from tadf_workbench.ui.angle_scan_table import AngleScanTable
from tadf_workbench.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {
        10: (3.4, 2.9), 30: (3.0, 2.95), 90: (3.4, 3.0),
    })))


def test_table_populates_rows(qtbot, scan):
    t = AngleScanTable()
    qtbot.addWidget(t)
    t.set_scan(scan)
    assert t.rowCount() == 3
    assert t.columnCount() == 8


def test_table_highlights_tadf_rows(qtbot, scan):
    t = AngleScanTable()
    qtbot.addWidget(t)
    t.set_scan(scan, threshold_ev=0.2)
    # 30° gap is 0.05, 10° gap is 0.5, 90° gap is 0.4
    assert t.is_row_highlighted(1)  # after sort: 10, 30, 90 → row 1 = 30°
    assert not t.is_row_highlighted(0)
    assert not t.is_row_highlighted(2)


def test_table_emits_angle_selected(qtbot, scan):
    t = AngleScanTable()
    qtbot.addWidget(t)
    t.set_scan(scan)
    received = []
    t.angle_selected.connect(received.append)
    t.selectRow(1)
    assert received and received[-1] == 30.0
