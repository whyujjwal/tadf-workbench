# tests/test_e2e_angle_scan.py
import pytest
pytest.importorskip("PyQt5")
from tadf_workbench.parsers import parse_angle_folder
from tadf_workbench.ui import AngleScanWindow
from tests.fixtures.build_fixtures import write_scan_folder


def test_full_flow(qtbot, tmp_path):
    folder = write_scan_folder(tmp_path / "demo", {
        10: (3.40, 2.90),  # gap 0.50
        20: (3.30, 2.92),  # gap 0.38
        30: (3.20, 2.95),  # gap 0.25
        40: (3.10, 2.98),  # gap 0.12 — TADF
        50: (3.05, 3.00),  # gap 0.05 — TADF (minimum)
        60: (3.08, 3.00),  # gap 0.08 — TADF
        70: (3.15, 2.99),  # gap 0.16 — TADF
        80: (3.25, 2.97),  # gap 0.28
        90: (3.40, 2.93),  # gap 0.47
    })
    scan = parse_angle_folder(str(folder))
    assert len(scan) == 9

    w = AngleScanWindow()
    qtbot.addWidget(w)
    w.set_scan(scan)

    assert w.table.rowCount() == 9
    # TADF candidates at θ = 40, 50, 60, 70
    assert sum(1 for r in range(w.table.rowCount()) if w.table.is_row_highlighted(r)) == 4

    # Clicking on θ=50 via the viewer updates status bar
    w.viewer_3d.set_angle(50)
    qtbot.wait(100)
    assert "50" in w.statusBar().currentMessage() or \
           "min gap" in w.statusBar().currentMessage().lower()
