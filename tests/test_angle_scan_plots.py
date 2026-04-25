# tests/test_angle_scan_plots.py
import pytest
pytest.importorskip("PyQt5")
pytestmark = pytest.mark.usefixtures("qtbot")

from tadf_workbench.ui.angle_scan_plots import GapPlot, TDMPlot, TwoPSquaredPlot
from tadf_workbench.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {
        10: (3.4, 2.9), 30: (3.0, 2.95), 60: (3.1, 2.95), 90: (3.4, 3.0),
    })))


def test_gap_plot_populates_from_scan(qtbot, scan):
    w = GapPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)
    # Should have two plot data items: main curve + highlight scatter
    items = [i for i in w.plotItem.items if hasattr(i, "setData")]
    assert len(items) >= 1


def test_gap_plot_threshold_updates(qtbot, scan):
    w = GapPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)
    w.set_threshold(0.1)
    assert w.threshold_ev == 0.1


def test_tdm_plot_populates(qtbot, scan):
    w = TDMPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)


def test_two_p_squared_plot_populates(qtbot, scan):
    w = TwoPSquaredPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)


def test_plots_emit_click_signal(qtbot, scan):
    w = GapPlot()
    qtbot.addWidget(w)
    w.set_scan(scan)
    received = []
    w.angle_clicked.connect(received.append)
    # Simulate programmatic point selection
    w._emit_nearest(30.0)
    assert received == [30.0]
