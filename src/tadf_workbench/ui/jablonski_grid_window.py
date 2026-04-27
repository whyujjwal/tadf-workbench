"""Grid window: one mini Jablonski diagram per θ, side by side.

Use this for spotting how the singlet/triplet manifold deforms across the
scan. Tiles are click-to-select — clicking any tile drives the main window
to that angle.
"""
from typing import Optional

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QFrame, QGridLayout, QHBoxLayout, QLabel, QMainWindow, QPushButton,
    QScrollArea, QVBoxLayout, QWidget,
)

from ..core import AngleScan
from .jablonski_panel import JablonskiPanel


_BG = "#0f1216"
_CARD = "#1a1e26"
_BORDER = "#2a313b"
_FG = "#eaf2fb"
_DIM = "#8a97a8"


class JablonskiGridWindow(QMainWindow):
    """Non-modal window showing every θ's Jablonski in a scrollable grid."""

    angle_picked = pyqtSignal(float)

    def __init__(self, scan: AngleScan, parent=None):
        super().__init__(parent)
        self.setWindowTitle(f"Jablonski Grid — {scan.name}")
        self.resize(1280, 820)
        self._scan = scan
        self._panels: list[JablonskiPanel] = []
        self._mode = "key"
        self.setStyleSheet(f"""
            QMainWindow, QWidget {{ background-color: {_BG}; color: {_FG}; }}
            QScrollArea {{ border: none; }}
            QPushButton {{
                background: {_CARD}; color: {_FG};
                border: 1px solid {_BORDER}; border-radius: 4px;
                padding: 4px 10px;
            }}
            QPushButton:checked {{
                background: #094771; border-color: #4fc3f7;
            }}
        """)
        self._build_ui()
        self._populate()

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # Header bar with mode toggle
        bar = QHBoxLayout()
        bar.setContentsMargins(2, 0, 2, 0)
        title = QLabel(
            f"Jablonski diagrams — {len(self._scan.points)} angles"
        )
        f = title.font(); f.setBold(True); f.setPointSize(11); title.setFont(f)
        bar.addWidget(title)
        bar.addStretch(1)

        mode_label = QLabel("Tile mode:")
        mode_label.setStyleSheet(f"color:{_DIM};")
        self.btn_key = QPushButton("Key only")
        self.btn_all = QPushButton("All states")
        self.btn_key.setCheckable(True); self.btn_all.setCheckable(True)
        self.btn_key.setChecked(True)
        self.btn_key.clicked.connect(lambda: self._set_mode("key"))
        self.btn_all.clicked.connect(lambda: self._set_mode("all"))
        bar.addWidget(mode_label)
        bar.addWidget(self.btn_key)
        bar.addWidget(self.btn_all)
        root.addLayout(bar)

        # Scrollable grid
        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        self._grid_host = QWidget()
        self._grid = QGridLayout(self._grid_host)
        self._grid.setContentsMargins(4, 4, 4, 4)
        self._grid.setSpacing(8)
        self._scroll.setWidget(self._grid_host)
        root.addWidget(self._scroll, 1)

    def _populate(self):
        # Clear previous (in case scan is re-set later)
        while self._grid.count():
            item = self._grid.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
        self._panels.clear()

        points = self._scan.sorted_points
        if not points:
            empty = QLabel("No angle points in this scan.")
            empty.setStyleSheet(f"color:{_DIM};")
            empty.setAlignment(Qt.AlignCenter)
            self._grid.addWidget(empty, 0, 0)
            return

        cols = 3
        for idx, point in enumerate(points):
            tile = self._make_tile(point)
            r, c = divmod(idx, cols)
            self._grid.addWidget(tile, r, c)

    def _make_tile(self, point) -> QWidget:
        card = QFrame()
        card.setStyleSheet(
            f"QFrame {{ background:{_CARD}; border:1px solid {_BORDER}; "
            f"border-radius:6px; }}"
        )
        v = QVBoxLayout(card)
        v.setContentsMargins(6, 6, 6, 6)
        v.setSpacing(4)

        panel = JablonskiPanel(show_mode_selector=False, compact=True)
        panel.setMinimumHeight(320)
        panel.set_mode(self._mode)
        panel.set_point(point)
        panel.angle_picked.connect(self.angle_picked.emit)
        self._panels.append(panel)
        v.addWidget(panel, 1)

        return card

    def _set_mode(self, mode: str):
        self._mode = mode
        self.btn_key.setChecked(mode == "key")
        self.btn_all.setChecked(mode == "all")
        for p in self._panels:
            p.set_mode(mode)
