"""2-D plot widgets for the angle-scan workspace."""
from typing import Optional

import numpy as np
import pyqtgraph as pg
from PyQt5.QtCore import pyqtSignal, Qt
from PyQt5.QtGui import QPen, QColor

from ..core import AngleScan
from ..analysis import gap_curve, tdm_curve, two_p_squared_curve


# Match the main-window dark theme
_BG = "#111518"
_FG = "#ddeeff"
_GRID = "#252d38"
_CURVE = "#4fc3f7"
_HIGHLIGHT = "#4caf50"
_THRESHOLD = "#ff8a65"


def _configure_axes(plot: pg.PlotWidget, *, xlabel: str, ylabel: str, title: str):
    plot.setBackground(_BG)
    plot.showGrid(x=True, y=True, alpha=0.15)
    plot.setLabel("bottom", xlabel, color=_FG)
    plot.setLabel("left", ylabel, color=_FG)
    plot.setTitle(title, color=_FG, size="11pt")
    for axis in ("bottom", "left"):
        ax = plot.getAxis(axis)
        ax.setPen(QColor(_GRID))
        ax.setTextPen(QColor(_FG))


class _AnglePlotBase(pg.PlotWidget):
    """Shared behaviour: set_scan, click-to-emit nearest angle."""

    angle_clicked = pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._scan: Optional[AngleScan] = None
        self._curve: Optional[pg.PlotDataItem] = None
        self.scene().sigMouseClicked.connect(self._on_scene_click)

    def set_scan(self, scan: AngleScan):
        self._scan = scan
        self._redraw()

    def _redraw(self):  # overridden
        raise NotImplementedError

    def _on_scene_click(self, ev):
        if not self._scan or not self._scan.points:
            return
        vb = self.plotItem.vb
        mouse_point = vb.mapSceneToView(ev.scenePos())
        nearest = self._scan.nearest_point(mouse_point.x())
        if nearest is not None:
            self.angle_clicked.emit(nearest.angle_deg)

    def _emit_nearest(self, angle_deg: float):
        """Test hook."""
        if self._scan is None:
            return
        nearest = self._scan.nearest_point(angle_deg)
        if nearest is not None:
            self.angle_clicked.emit(nearest.angle_deg)


class GapPlot(_AnglePlotBase):
    def __init__(self, parent=None):
        super().__init__(parent)
        _configure_axes(self, xlabel="θ (deg)", ylabel="ΔE(S₁−T₁) (eV)",
                        title="Singlet–Triplet Gap")
        self.threshold_ev = 0.2
        self._threshold_line: Optional[pg.InfiniteLine] = None
        self._highlight: Optional[pg.ScatterPlotItem] = None

    def set_threshold(self, value: float):
        self.threshold_ev = value
        self._redraw()

    def _redraw(self):
        self.clear()
        self._threshold_line = None
        self._highlight = None
        if not self._scan or not self._scan.points:
            return
        xs, ys = gap_curve(self._scan)
        pen = pg.mkPen(_CURVE, width=2)
        self._curve = self.plot(xs, ys, pen=pen, symbol="o",
                                symbolBrush=_CURVE, symbolSize=8)
        # Dashed threshold line
        thr_pen = pg.mkPen(_THRESHOLD, width=1.5, style=Qt.DashLine)
        self._threshold_line = pg.InfiniteLine(
            pos=self.threshold_ev, angle=0, pen=thr_pen,
            label=f"{self.threshold_ev:.2f} eV",
            labelOpts={"color": _THRESHOLD, "position": 0.95},
        )
        self.addItem(self._threshold_line)
        # Highlight TADF candidates (gap ≤ threshold)
        xs_h = [x for x, y in zip(xs, ys) if y is not None and y <= self.threshold_ev]
        ys_h = [y for y in ys if y is not None and y <= self.threshold_ev]
        if xs_h:
            self._highlight = pg.ScatterPlotItem(
                x=xs_h, y=ys_h,
                brush=pg.mkBrush(_HIGHLIGHT), pen=pg.mkPen("w", width=1.5),
                size=14, symbol="o",
            )
            self.addItem(self._highlight)


class TDMPlot(_AnglePlotBase):
    def __init__(self, parent=None):
        super().__init__(parent)
        _configure_axes(self, xlabel="θ (deg)", ylabel="|TDM(S₁)| (Au)",
                        title="S₁ Transition Dipole Moment")

    def _redraw(self):
        self.clear()
        if not self._scan or not self._scan.points:
            return
        xs, ys = tdm_curve(self._scan)
        self._curve = self.plot(xs, ys, pen=pg.mkPen("#ffeb3b", width=2),
                                symbol="s", symbolBrush="#ffeb3b", symbolSize=8)


class TwoPSquaredPlot(_AnglePlotBase):
    def __init__(self, parent=None):
        super().__init__(parent)
        _configure_axes(self, xlabel="θ (deg)", ylabel="2·P²",
                        title="Dominant S₁ Configuration Weight (2·P²)")

    def _redraw(self):
        self.clear()
        if not self._scan or not self._scan.points:
            return
        xs, ys = two_p_squared_curve(self._scan)
        self._curve = self.plot(xs, ys, pen=pg.mkPen("#b388ff", width=2),
                                symbol="t", symbolBrush="#b388ff", symbolSize=9)
