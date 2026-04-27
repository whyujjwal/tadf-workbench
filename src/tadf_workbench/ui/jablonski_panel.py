"""Jablonski energy-level diagram for one selected angle.

Two columns (singlets left, triplets right) on a shared eV y-axis. Highlights
S₁, the f-max singlet (Sₙ), and T₁, with multi-line labels showing the
dominant orbital transition and 2·P² as a percentage. Toggle between a
key-states-only and an all-states view.
"""
from typing import List, Optional

import pyqtgraph as pg
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QFont, QPen
from PyQt5.QtWidgets import (
    QButtonGroup, QHBoxLayout, QLabel, QRadioButton, QSizePolicy,
    QVBoxLayout, QWidget,
)

from ..core import AnglePoint
from ..analysis import JablonskiLayout, jablonski_layout, empty_layout


# Match the workbench dark palette
_BG = "#111518"
_FG = "#eaf2fb"
_DIM = "#8a97a8"
_GRID = "#252d38"

_SINGLET = "#4fc3f7"
_SINGLET_BRIGHT = "#7dd3fc"   # for the f-max singlet (Sn) — slightly different
_TRIPLET = "#ff8a65"
_GROUND = "#cfd8dc"

_FAINT = "#39424f"             # non-highlighted levels in "all" mode
_ABS_ARROW = "#ffeb3b"
_ISC_ARROW = "#b388ff"

_X_SINGLET = 1.5
_X_TRIPLET = 3.5
_LEVEL_HALF_WIDTH = 0.55
_X_RANGE = (-1.4, 6.2)  # extra room on each side for labels


class JablonskiPanel(QWidget):
    """Live single-angle Jablonski. set_point(point) updates the drawing."""

    angle_picked = pyqtSignal(float)  # for click-to-select inside grid view

    def __init__(self, parent=None, *, show_mode_selector: bool = True,
                 compact: bool = False):
        super().__init__(parent)
        self._point: Optional[AnglePoint] = None
        self._mode: str = "key"
        self._compact = compact
        self._show_mode_selector = show_mode_selector
        self._build_ui()

    # ── UI ────────────────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(2, 2, 2, 2)
        root.setSpacing(4)

        if self._show_mode_selector and not self._compact:
            mode_row = QHBoxLayout()
            mode_row.setContentsMargins(8, 4, 8, 0)
            label = QLabel("States:")
            label.setStyleSheet(f"color:{_DIM}; font-size:11px;")
            self.rb_key = QRadioButton("Key only (S₀, S₁, Sₙ, T₁)")
            self.rb_all = QRadioButton("All states")
            self.rb_key.setChecked(True)
            for rb in (self.rb_key, self.rb_all):
                rb.setStyleSheet(f"color:{_FG}; font-size:11px; padding:0 6px;")
            self._mode_group = QButtonGroup(self)
            self._mode_group.addButton(self.rb_key, 0)
            self._mode_group.addButton(self.rb_all, 1)
            self._mode_group.buttonClicked.connect(self._on_mode_changed)
            mode_row.addWidget(label)
            mode_row.addWidget(self.rb_key)
            mode_row.addWidget(self.rb_all)
            mode_row.addStretch(1)
            root.addLayout(mode_row)

        self.title_label = QLabel("Jablonski — no scan loaded")
        self.title_label.setAlignment(Qt.AlignCenter)
        title_font = QFont()
        title_font.setBold(True)
        title_font.setPointSize(10 if self._compact else 11)
        self.title_label.setFont(title_font)
        self.title_label.setStyleSheet(
            f"color:{_FG}; padding:{2 if self._compact else 4}px;")
        root.addWidget(self.title_label)

        self.plot = pg.PlotWidget()
        self.plot.setBackground(_BG)
        self.plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.plot.showGrid(x=False, y=True, alpha=0.10)
        self.plot.setMouseEnabled(x=False, y=False)
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        # X axis is categorical — hide
        self.plot.getAxis("bottom").setStyle(showValues=False)
        self.plot.getAxis("bottom").setPen(QColor(_GRID))
        self.plot.getAxis("left").setPen(QColor(_GRID))
        self.plot.getAxis("left").setTextPen(QColor(_FG))
        self.plot.setLabel("left", "Energy (eV)", color=_FG)
        self.plot.setXRange(_X_RANGE[0], _X_RANGE[1], padding=0)
        self.plot.scene().sigMouseClicked.connect(self._on_scene_click)
        root.addWidget(self.plot, 1)

    # ── Public API ────────────────────────────────────────────────────

    def set_point(self, point: Optional[AnglePoint]):
        self._point = point
        self._redraw()

    def set_mode(self, mode: str):
        if mode not in ("key", "all"):
            return
        self._mode = mode
        if hasattr(self, "rb_key") and hasattr(self, "rb_all"):
            self.rb_key.setChecked(mode == "key")
            self.rb_all.setChecked(mode == "all")
        self._redraw()

    # ── Handlers ─────────────────────────────────────────────────────

    def _on_mode_changed(self, _btn):
        self._mode = "key" if self.rb_key.isChecked() else "all"
        self._redraw()

    def _on_scene_click(self, _ev):
        if self._point is not None:
            self.angle_picked.emit(self._point.angle_deg)

    # ── Rendering ────────────────────────────────────────────────────

    def _redraw(self):
        self.plot.clear()
        if self._point is None:
            self.title_label.setText("Jablonski — no point selected")
            return

        layout = jablonski_layout(self._point, mode=self._mode)
        self.title_label.setText(layout.title)

        if not layout.has_data:
            return

        max_e = max(layout.max_energy_ev, 1.0)
        self.plot.setYRange(-0.15, max_e * 1.15, padding=0)

        # Column header text
        self._add_column_header("Singlets", _X_SINGLET, max_e * 1.10)
        self._add_column_header("Triplets", _X_TRIPLET, max_e * 1.10)

        # Levels (faint first so highlighted draw on top)
        ordered = sorted(layout.levels, key=lambda lv: lv.is_highlighted)
        for lv in ordered:
            self._draw_level(lv)

        # Arrows on top
        for ar in layout.arrows:
            self._draw_arrow(ar)

    def _add_column_header(self, text: str, x: float, y: float):
        item = pg.TextItem(text, color=_DIM, anchor=(0.5, 0.5))
        f = QFont(); f.setBold(True); f.setPointSize(9)
        item.setFont(f)
        item.setPos(x, y)
        self.plot.addItem(item)

    def _draw_level(self, lv):
        x = _X_SINGLET if lv.column == "singlet" else _X_TRIPLET
        x0 = x - _LEVEL_HALF_WIDTH
        x1 = x + _LEVEL_HALF_WIDTH

        if lv.role == "S0":
            color = _GROUND
            width = 3
        elif lv.role == "Sn":
            color = _SINGLET_BRIGHT
            width = 4
        elif lv.is_highlighted and lv.column == "singlet":
            color = _SINGLET
            width = 4
        elif lv.is_highlighted and lv.column == "triplet":
            color = _TRIPLET
            width = 4
        else:
            color = _FAINT
            width = 1

        pen = pg.mkPen(color, width=width)
        self.plot.plot([x0, x1], [lv.energy_ev, lv.energy_ev],
                       pen=pen, antialias=True)

        if lv.label:
            # Place label on the outer side of the column (to avoid the
            # absorption arrow that sits between S₀ and the level)
            if lv.column == "singlet":
                label_x = x0 - 0.05
                anchor = (1.0, 0.5)
            else:
                label_x = x1 + 0.05
                anchor = (0.0, 0.5)
            text_color = color if lv.is_highlighted else _DIM
            item = pg.TextItem(lv.label, color=text_color, anchor=anchor)
            f = QFont(); f.setPointSize(8 if self._compact else 9)
            item.setFont(f)
            item.setPos(label_x, lv.energy_ev)
            self.plot.addItem(item)

    def _draw_arrow(self, ar):
        if ar.kind == "absorption":
            x = _X_SINGLET
            color = _ABS_ARROW
        elif ar.kind == "isc":
            # Horizontal-ish arrow between S₁ on left and T₁ on right.
            self._draw_isc_arrow(ar)
            return
        else:
            return

        # Draw a vertical line and arrowhead at the top
        line_pen = pg.mkPen(color, width=2)
        self.plot.plot([x, x], [ar.from_energy_ev, ar.to_energy_ev],
                       pen=line_pen, antialias=True)
        head = pg.ArrowItem(
            pos=(x, ar.to_energy_ev), angle=-90,
            tipAngle=35, headLen=14, brush=color, pen=line_pen,
        )
        self.plot.addItem(head)
        # Label sits half-way along the arrow, slightly offset from the line
        if ar.label:
            mid_y = (ar.from_energy_ev + ar.to_energy_ev) / 2.0
            txt = pg.TextItem(ar.label, color=color, anchor=(0.0, 0.5))
            f = QFont(); f.setPointSize(8 if self._compact else 9)
            txt.setFont(f)
            txt.setPos(x + 0.06, mid_y)
            self.plot.addItem(txt)

    def _draw_isc_arrow(self, ar):
        # Wavy-ish horizontal arrow between S₁ and T₁ at their respective
        # energies.  Drawn as a single straight segment for clarity.
        color = _ISC_ARROW
        pen = pg.mkPen(color, width=2, style=Qt.DashLine)
        x1, x2 = _X_SINGLET + _LEVEL_HALF_WIDTH, _X_TRIPLET - _LEVEL_HALF_WIDTH
        self.plot.plot([x1, x2], [ar.from_energy_ev, ar.to_energy_ev],
                       pen=pen, antialias=True)
        # Arrow head pointing toward T₁ side
        # Compute angle so the head matches the line direction.
        dx = x2 - x1
        dy = ar.to_energy_ev - ar.from_energy_ev
        # pyqtgraph arrow's `angle` is the direction the tail points (deg).
        # We want the head to point from S₁ toward T₁.
        import math
        angle_deg = math.degrees(math.atan2(dy, dx))
        head = pg.ArrowItem(
            pos=(x2, ar.to_energy_ev), angle=180 - angle_deg,
            tipAngle=30, headLen=12, brush=color, pen=pen,
        )
        self.plot.addItem(head)

        if ar.label:
            mx = (x1 + x2) / 2.0
            my = (ar.from_energy_ev + ar.to_energy_ev) / 2.0
            txt = pg.TextItem(ar.label, color=color, anchor=(0.5, 1.1))
            f = QFont(); f.setPointSize(8 if self._compact else 9)
            txt.setFont(f)
            txt.setPos(mx, my)
            self.plot.addItem(txt)
