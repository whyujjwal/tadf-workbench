"""Tabular view of AngleScan metrics with TADF-row highlighting."""
from typing import Optional

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QBrush, QFont
from PyQt5.QtWidgets import QTableWidget, QTableWidgetItem, QAbstractItemView

from ..core import AngleScan
from ..analysis import summary_row, percent_2p2


_TADF_BG = QColor(40, 95, 45)
_TADF_FG = QColor("#c8ffd4")
_NORMAL_BG = QColor("#1e1e1e")
_NORMAL_FG = QColor("#d4d4d4")

_COLUMNS = [
    ("θ (deg)", "angle_deg"),
    ("E(S₁) eV", "s1_energy_ev"),
    ("E(T₁) eV", "t1_energy_ev"),
    ("ΔE eV", "gap_ev"),
    ("|TDM| Au", "tdm_magnitude"),
    ("2·P²", "two_p_squared"),
    ("S₁ dominant (i→a, %)", "s1_dominant"),
    ("S₁ 2nd jump (i→a, %)", "s1_second"),
]


def _fmt(val, key):
    if val is None:
        return "—"
    if key in ("s1_dominant", "s1_second"):
        i, a, c = val
        return f"{i} → {a}  ({percent_2p2(c):.1f}%)"
    if key == "angle_deg":
        return f"{val:.1f}"
    if isinstance(val, float):
        return f"{val:.4f}"
    return str(val)


class AngleScanTable(QTableWidget):
    angle_selected = pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._scan: Optional[AngleScan] = None
        self._threshold = 0.2
        self._highlight_rows: set = set()
        self.setColumnCount(len(_COLUMNS))
        self.setHorizontalHeaderLabels([c[0] for c in _COLUMNS])
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.setAlternatingRowColors(True)
        self.setStyleSheet("""
            QTableWidget {
                background-color: #1e1e1e; color: #d4d4d4;
                gridline-color: #333; alternate-background-color: #242424;
                selection-background-color: #094771;
            }
            QHeaderView::section {
                background-color: #2d2d30; color: #ddd;
                padding: 4px; border: 1px solid #444;
            }
        """)
        self.itemSelectionChanged.connect(self._emit_selected)

    def set_scan(self, scan: AngleScan, threshold_ev: float = 0.2):
        self._scan = scan
        self._threshold = threshold_ev
        self._highlight_rows.clear()
        points = scan.sorted_points
        self.setRowCount(len(points))
        for row, p in enumerate(points):
            row_data = summary_row(p)
            is_tadf = (p.s1_t1_gap_ev is not None and p.s1_t1_gap_ev <= threshold_ev)
            if is_tadf:
                self._highlight_rows.add(row)
            for col, (_, key) in enumerate(_COLUMNS):
                item = QTableWidgetItem(_fmt(row_data.get(key), key))
                item.setTextAlignment(Qt.AlignCenter)
                if is_tadf:
                    item.setBackground(QBrush(_TADF_BG))
                    item.setForeground(QBrush(_TADF_FG))
                    f = item.font(); f.setBold(True); item.setFont(f)
                self.setItem(row, col, item)
        self.resizeColumnsToContents()

    def is_row_highlighted(self, row: int) -> bool:
        return row in self._highlight_rows

    def _emit_selected(self):
        if not self._scan:
            return
        rows = self.selectionModel().selectedRows()
        if not rows:
            return
        idx = rows[0].row()
        points = self._scan.sorted_points
        if 0 <= idx < len(points):
            self.angle_selected.emit(points[idx].angle_deg)
