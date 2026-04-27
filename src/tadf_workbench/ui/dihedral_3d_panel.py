"""3D viewer for the current angle point.

Renders the molecule with two visually distinct fragments (donor / acceptor),
the rotatable bond highlighted, both fragment best-fit planes drawn as
translucent quads, and the IUPAC dihedral angle as an arc + degree label.
"""
from typing import Optional

import numpy as np
import pyqtgraph.opengl as gl
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QFont, QVector3D
from PyQt5.QtWidgets import (
    QCheckBox, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QSlider, QFrame,
    QSizePolicy,
)

from ..core import AngleScan, AnglePoint, Stage
from ..core.fragments import dihedral_angle, fit_plane, plane_quad
from ..utils import get_atom_color, get_atom_radius


# ── Colour palette ──
DONOR_TINT = np.array([0.40, 0.95, 0.50])   # green
ACCEPTOR_TINT = np.array([0.95, 0.55, 0.35])  # orange
ROTATABLE_BOND_COLOR = (1.0, 0.95, 0.30, 1.0)  # bright yellow
NEUTRAL_BOND_COLOR = (0.78, 0.82, 0.92, 0.95)


def _blend(element_rgb: tuple, tint: np.ndarray, weight: float = 0.55) -> tuple:
    """Blend an RGB tuple toward a tint colour, returning (r, g, b, a)."""
    base = np.array(element_rgb[:3])
    blended = base * (1 - weight) + tint * weight
    return (float(blended[0]), float(blended[1]), float(blended[2]), 1.0)


class Dihedral3DPanel(QWidget):
    angle_changed = pyqtSignal(float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._scan: Optional[AngleScan] = None
        self._current: Optional[AnglePoint] = None
        self._gl_items: list = []
        self._camera_fitted = False
        self._show_planes = True
        self._show_arc = True
        self._show_atom_numbers = False
        self._build_ui()

    def _build_ui(self):
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)

        self.view = gl.GLViewWidget()
        self.view.setBackgroundColor("#0d1014")
        self.view.opts["distance"] = 25
        self.view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        lay.addWidget(self.view, 1)

        # Bottom bar: slider + dihedral readout + toggles
        bar = QFrame()
        bar.setFixedHeight(86)
        bar.setStyleSheet("""
            QFrame { background: qlineargradient(x1:0,y1:0,x2:0,y2:1,
                      stop:0 #1a1e26, stop:1 #141820);
                     border-top: 1px solid #2f3640; }
            QLabel { background: transparent; }
        """)
        bar_lay = QVBoxLayout(bar)
        bar_lay.setContentsMargins(14, 6, 14, 6)
        bar_lay.setSpacing(2)

        # Top row: scan θ + measured dihedral + scan info
        top = QHBoxLayout()
        top.setSpacing(14)
        theta_tag = QLabel("SCAN θ")
        theta_tag.setStyleSheet("color:#8a97a8; font-size:10px; "
                                "font-weight:600; letter-spacing:1px;")
        top.addWidget(theta_tag)
        self.angle_label = QLabel("—°")
        self.angle_label.setStyleSheet("color:#4fc3f7; font-size:18px; "
                                       "font-weight:700;")
        top.addWidget(self.angle_label)

        sep = QLabel("│")
        sep.setStyleSheet("color:#2f3640;")
        top.addWidget(sep)

        meas_tag = QLabel("MEASURED N₁–B₁–B₂–N₂")
        meas_tag.setStyleSheet("color:#8a97a8; font-size:10px; "
                               "font-weight:600; letter-spacing:1px;")
        top.addWidget(meas_tag)
        self.dihedral_label = QLabel("—°")
        self.dihedral_label.setStyleSheet("color:#ffeb3b; font-size:14px; "
                                          "font-weight:700;")
        top.addWidget(self.dihedral_label)

        top.addStretch()
        self.scan_info_label = QLabel("")
        self.scan_info_label.setStyleSheet("color:#8a97a8; font-size:10px;")
        top.addWidget(self.scan_info_label)
        bar_lay.addLayout(top)

        # Slider
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(180)
        self.slider.setTickPosition(QSlider.TicksBelow)
        self.slider.setTickInterval(10)
        self.slider.setStyleSheet("""
            QSlider::groove:horizontal {
                height: 6px; background: #2a313b;
                border-radius: 3px; margin: 0 2px;
            }
            QSlider::sub-page:horizontal {
                background: qlineargradient(x1:0,y1:0,x2:1,y2:0,
                    stop:0 #4fc3f7, stop:1 #1177bb);
                border-radius: 3px;
            }
            QSlider::handle:horizontal {
                background: #4fc3f7; border: 2px solid #eaf2fb;
                width: 16px; height: 16px; margin: -6px 0;
                border-radius: 9px;
            }
            QSlider::handle:horizontal:hover {
                background: #81d4fa; border: 2px solid white;
            }
        """)
        self.slider.valueChanged.connect(self._on_slider_change)
        bar_lay.addWidget(self.slider)

        # Legend row
        legend = QHBoxLayout()
        legend.setSpacing(12)
        legend.addWidget(self._make_swatch("Donor fragment", "#66bb6a"))
        legend.addWidget(self._make_swatch("Acceptor fragment", "#ff8a65"))
        legend.addWidget(self._make_swatch("Rotatable bond", "#ffeb3b"))
        legend.addStretch()
        self._numbers_chk = QCheckBox("Atom #")
        self._numbers_chk.setStyleSheet(
            "QCheckBox { color:#b5c2d1; font-size:10px; }"
            "QCheckBox::indicator { width:12px; height:12px; }"
        )
        self._numbers_chk.toggled.connect(self._on_numbers_toggled)
        legend.addWidget(self._numbers_chk)
        bar_lay.addLayout(legend)

        lay.addWidget(bar)

    def _make_swatch(self, text: str, color: str) -> QWidget:
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.setSpacing(5)
        dot = QLabel("●")
        dot.setStyleSheet(f"color:{color}; font-size:14px;")
        lab = QLabel(text)
        lab.setStyleSheet("color:#b5c2d1; font-size:10px;")
        h.addWidget(dot)
        h.addWidget(lab)
        return w

    # ── Public API ──

    @property
    def current_angle(self) -> Optional[float]:
        return self._current.angle_deg if self._current else None

    def set_scan(self, scan: AngleScan):
        self._scan = scan
        self._camera_fitted = False
        if not scan.points:
            return
        pts = scan.sorted_points
        lo, hi = int(pts[0].angle_deg), int(pts[-1].angle_deg)
        self.slider.blockSignals(True)
        self.slider.setMinimum(lo)
        self.slider.setMaximum(hi)
        self.slider.blockSignals(False)
        self.scan_info_label.setText(
            f"{len(pts)} points · {lo}°–{hi}° · "
            f"atoms {pts[0].molecule.current_stage.atom_count}"
        )
        self.set_angle(pts[0].angle_deg)

    def set_angle(self, angle_deg: float):
        if not self._scan or not self._scan.points:
            return
        target = self._scan.nearest_point(angle_deg)
        if target is None:
            return
        if self._current is not None and target.angle_deg == self._current.angle_deg:
            return
        self._current = target
        self.slider.blockSignals(True)
        self.slider.setValue(int(round(target.angle_deg)))
        self.slider.blockSignals(False)
        self.angle_label.setText(f"{target.angle_deg:.1f}°")
        self._render_point(target)
        self.angle_changed.emit(target.angle_deg)

    # ── Internals ──

    def _on_slider_change(self, value: int):
        self.set_angle(float(value))

    def _on_numbers_toggled(self, checked: bool):
        self._show_atom_numbers = checked
        if self._current is not None:
            self._render_point(self._current)

    def _clear_items(self):
        for it in self._gl_items:
            self.view.removeItem(it)
        self._gl_items.clear()

    def _render_point(self, point: AnglePoint):
        self._clear_items()
        stage = point.molecule.current_stage
        if stage is None or stage.atom_count == 0:
            return

        donor = self._scan.donor_fragment if self._scan else None
        accept = self._scan.acceptor_fragment if self._scan else None
        bond = self._scan.rotatable_bond if self._scan else None

        if not self._camera_fitted:
            self._fit_camera_to_stage(stage)
            self._camera_fitted = True

        self._draw_atoms(stage, donor, accept)
        self._draw_bonds(stage, bond)
        if self._show_atom_numbers:
            self._draw_atom_numbers(stage)

        if bond is not None and donor is not None and accept is not None:
            self._draw_rotatable_axis(stage, bond)
            if self._show_planes:
                self._draw_fragment_plane(stage, donor, bond, DONOR_TINT)
                self._draw_fragment_plane(stage, accept, bond, ACCEPTOR_TINT)
            self._draw_dihedral(stage)
        else:
            self.dihedral_label.setText("(set rotatable bond)")

    def _fit_camera_to_stage(self, stage: Stage):
        positions = stage.positions
        bbox_min = positions.min(axis=0)
        bbox_max = positions.max(axis=0)
        center = (bbox_min + bbox_max) / 2
        extent = float(np.linalg.norm(bbox_max - bbox_min))
        distance = max(12.0, extent * 1.8)
        self.view.opts["center"] = QVector3D(*center)
        self.view.setCameraPosition(distance=distance, elevation=18, azimuth=30)

    def _draw_atoms(self, stage: Stage, donor: Optional[set],
                    accept: Optional[set]):
        md = gl.MeshData.sphere(rows=16, cols=28)
        for atom in stage.atoms:
            base = get_atom_color(atom.atomic_number)
            base = base if len(base) == 4 else (*base, 1.0)
            if donor is not None and atom.index in donor:
                color = _blend(base, DONOR_TINT, weight=0.45)
            elif accept is not None and atom.index in accept:
                color = _blend(base, ACCEPTOR_TINT, weight=0.45)
            else:
                color = base
            radius = get_atom_radius(atom.atomic_number) * 0.55
            item = gl.GLMeshItem(meshdata=md, smooth=True, color=color,
                                 shader="shaded", glOptions="opaque")
            item.scale(radius, radius, radius)
            item.translate(atom.x, atom.y, atom.z)
            self.view.addItem(item)
            self._gl_items.append(item)

    def _draw_atom_numbers(self, stage: Stage):
        """Render small index labels next to each atom — toggled via checkbox.

        Helps the user see which atom corresponds to which index when picking
        the rotatable bond manually.
        """
        font = QFont()
        font.setPointSize(9)
        font.setBold(True)
        for atom in stage.atoms:
            radius = get_atom_radius(atom.atomic_number) * 0.55
            offset = radius + 0.12
            try:
                txt = gl.GLTextItem(
                    pos=(atom.x + offset, atom.y + offset, atom.z + offset),
                    text=str(atom.index),
                    color=QColor("#ffd54f"),
                    font=font,
                )
            except TypeError:
                # Older pyqtgraph: GLTextItem doesn't take font kwarg
                txt = gl.GLTextItem(
                    pos=(atom.x + offset, atom.y + offset, atom.z + offset),
                    text=str(atom.index),
                    color=QColor("#ffd54f"),
                )
            self.view.addItem(txt)
            self._gl_items.append(txt)

    def _draw_bonds(self, stage: Stage, rotatable_bond: Optional[tuple]):
        bonds = stage.find_bonds()
        if not bonds:
            return
        positions = stage.positions
        # Render the rotatable bond separately (thicker / coloured)
        normal_pts = []
        for i, j in bonds:
            if rotatable_bond is not None and {i, j} == set(rotatable_bond):
                continue
            normal_pts.append(positions[i])
            normal_pts.append(positions[j])
        if normal_pts:
            line = gl.GLLinePlotItem(
                pos=np.array(normal_pts), color=NEUTRAL_BOND_COLOR,
                width=3.0, mode="lines", antialias=True,
            )
            self.view.addItem(line)
            self._gl_items.append(line)

    def _draw_rotatable_axis(self, stage: Stage, bond: tuple):
        positions = stage.positions
        b1, b2 = bond
        p1, p2 = positions[b1], positions[b2]
        # Thick yellow segment between the two bond atoms
        seg = gl.GLLinePlotItem(
            pos=np.array([p1, p2]), color=ROTATABLE_BOND_COLOR,
            width=8.0, mode="lines", antialias=True,
        )
        self.view.addItem(seg)
        self._gl_items.append(seg)
        # Extend the rotation axis a bit beyond each atom for clarity
        axis = (p2 - p1)
        axis_norm = axis / (np.linalg.norm(axis) + 1e-9)
        ext = axis_norm * 1.6
        ghost = gl.GLLinePlotItem(
            pos=np.array([p1 - ext, p1, p2, p2 + ext]),
            color=(1.0, 0.95, 0.30, 0.35),
            width=2.0, mode="lines", antialias=True,
        )
        self.view.addItem(ghost)
        self._gl_items.append(ghost)

    def _draw_fragment_plane(self, stage: Stage, fragment: set,
                             bond: tuple, tint: np.ndarray):
        if len(fragment) < 3:
            return
        positions = stage.positions
        frag_pts = positions[list(fragment)]
        centroid, normal = fit_plane(frag_pts)
        # In-plane hint: the vector from bond midpoint toward the fragment centroid
        b1, b2 = bond
        bond_mid = (positions[b1] + positions[b2]) / 2
        hint = centroid - bond_mid
        # Size the quad to roughly span the fragment
        spread = float(np.linalg.norm(frag_pts.std(axis=0))) * 2.2 + 1.5
        quad = plane_quad(centroid, normal, half_extent=spread, in_plane_hint=hint)
        # Draw as a translucent surface (mesh of two triangles)
        verts = quad
        faces = np.array([[0, 1, 2], [0, 2, 3]])
        color = (*tint.tolist(), 0.18)
        mesh_data = gl.MeshData(vertexes=verts, faces=faces)
        mesh = gl.GLMeshItem(
            meshdata=mesh_data, smooth=False, color=color,
            shader="balloon", glOptions="translucent", drawEdges=True,
            edgeColor=(*tint.tolist(), 0.55),
        )
        self.view.addItem(mesh)
        self._gl_items.append(mesh)

    def _draw_dihedral(self, stage: Stage):
        scan = self._scan
        if scan is None or scan.rotatable_bond is None:
            return
        b1, b2 = scan.rotatable_bond
        n1 = scan.donor_reference_atom
        n2 = scan.acceptor_reference_atom
        if n1 is None or n2 is None:
            self.dihedral_label.setText("(no reference atoms)")
            return
        positions = stage.positions
        p_n1, p_b1, p_b2, p_n2 = positions[n1], positions[b1], positions[b2], positions[n2]
        angle = dihedral_angle(p_n1, p_b1, p_b2, p_n2)
        self.dihedral_label.setText(f"{angle:+.1f}°")

        if not self._show_arc:
            return

        # Build an arc in the plane perpendicular to the bond axis.
        bond_mid = (p_b1 + p_b2) / 2
        axis = (p_b2 - p_b1)
        axis = axis / (np.linalg.norm(axis) + 1e-9)
        # Project (n1 - b1) and (n2 - b2) onto the plane perpendicular to axis
        v1 = (p_n1 - p_b1) - np.dot(p_n1 - p_b1, axis) * axis
        v2 = (p_n2 - p_b2) - np.dot(p_n2 - p_b2, axis) * axis
        if np.linalg.norm(v1) < 1e-6 or np.linalg.norm(v2) < 1e-6:
            return
        v1u = v1 / np.linalg.norm(v1)
        v2u = v2 / np.linalg.norm(v2)
        # Place arc near the bond midpoint
        radius = max(0.6, min(np.linalg.norm(v1), np.linalg.norm(v2)) * 0.5)
        steps = 36
        # Unsigned angle [0, 180]
        cos_t = float(np.clip(np.dot(v1u, v2u), -1.0, 1.0))
        unsigned = np.arccos(cos_t)
        # Build orthonormal basis (v1u, w) in the perpendicular plane
        w = v2u - np.dot(v2u, v1u) * v1u
        if np.linalg.norm(w) < 1e-6:
            return
        w = w / np.linalg.norm(w)
        ts = np.linspace(0, unsigned, steps)
        arc = np.stack([
            bond_mid + radius * (np.cos(t) * v1u + np.sin(t) * w)
            for t in ts
        ])
        arc_item = gl.GLLinePlotItem(
            pos=arc, color=(1.0, 0.95, 0.30, 0.85),
            width=3.0, mode="line_strip", antialias=True,
        )
        self.view.addItem(arc_item)
        self._gl_items.append(arc_item)

        # Also draw the two reference vectors as thin guide lines
        guides = np.array([p_b1, p_n1, p_b2, p_n2])
        guide_item = gl.GLLinePlotItem(
            pos=guides, color=(1.0, 0.95, 0.30, 0.55),
            width=1.5, mode="lines", antialias=True,
        )
        self.view.addItem(guide_item)
        self._gl_items.append(guide_item)
