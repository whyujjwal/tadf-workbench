"""Top-level angle-scan analysis workspace — molecule-centric layout."""
from typing import Optional

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QKeySequence
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QSplitter, QAction,
    QFileDialog, QInputDialog, QMessageBox, QLabel, QFrame, QSizePolicy,
)

from ..core import AngleScan, AnglePoint
from ..parsers import parse_angle_folder, angle_from_filename
from ..parsers.gaussian_parser import parse_gaussian_file
from .angle_scan_plots import GapPlot, TDMPlot, TwoPSquaredPlot
from .angle_scan_table import AngleScanTable
from .dihedral_3d_panel import Dihedral3DPanel
from .donor_acceptor_dialog import DonorAcceptorDialog


# ── Palette ─────────────────────────────────────────────────────────────────
BG_WINDOW = "#0f1216"
BG_CARD = "#1a1e26"
BG_CARD_HI = "#232833"
BORDER = "#2a313b"
TEXT_DIM = "#8a97a8"
TEXT_MED = "#b5c2d1"
TEXT_BRIGHT = "#eaf2fb"
ACCENT = "#4fc3f7"
SINGLET = "#4fc3f7"
TRIPLET = "#ff8a65"
GOOD = "#66bb6a"
BAD = "#ef5350"
WARN = "#ffb74d"


# ── Tiny presentational widgets ──────────────────────────────────────────────

class StatChip(QFrame):
    """Horizontal chip: TITLE value subtitle — compact for the top bar."""

    def __init__(self, title: str, parent=None):
        super().__init__(parent)
        self.setStyleSheet(f"""
            StatChip {{
                background: {BG_CARD}; border: 1px solid {BORDER};
                border-radius: 6px;
            }}
        """)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(12, 6, 12, 6)
        lay.setSpacing(10)
        self._title = QLabel(title)
        self._title.setStyleSheet(f"color:{TEXT_DIM}; font-size:10px; "
                                  f"font-weight:600; letter-spacing:1px;")
        self._value = QLabel("—")
        self._value.setStyleSheet(f"color:{TEXT_BRIGHT}; font-size:16px; "
                                  f"font-weight:700;")
        self._sub = QLabel("")
        self._sub.setStyleSheet(f"color:{TEXT_DIM}; font-size:10px;")
        lay.addWidget(self._title)
        lay.addWidget(self._value)
        lay.addWidget(self._sub, 1)
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)

    def set_value(self, text: str, subtitle: str = "",
                  tint: Optional[str] = None):
        self._value.setText(text)
        self._sub.setText(subtitle)
        color = tint or TEXT_BRIGHT
        self._value.setStyleSheet(f"color:{color}; font-size:16px; "
                                  f"font-weight:700;")


class CurrentStateCard(QFrame):
    """Compact current-θ details — fits into a narrow sidebar column."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setStyleSheet(f"""
            CurrentStateCard {{
                background: {BG_CARD}; border: 1px solid {BORDER};
                border-radius: 6px;
            }}
            QLabel {{ background: transparent; }}
        """)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(12, 10, 12, 10)
        lay.setSpacing(6)

        # S1 / T1 side-by-side (compact)
        row = QHBoxLayout()
        row.setSpacing(8)
        self.s1_box = self._make_state_box("S₁", SINGLET)
        self.t1_box = self._make_state_box("T₁", TRIPLET)
        row.addWidget(self.s1_box["frame"])
        row.addWidget(self.t1_box["frame"])
        lay.addLayout(row)

        self.gap_label = QLabel("ΔE = —")
        self.gap_label.setAlignment(Qt.AlignCenter)
        self.gap_label.setStyleSheet(
            f"color:{TEXT_BRIGHT}; font-size:13px; font-weight:700; "
            f"padding:6px; border-radius:4px; background:{BG_CARD_HI};"
        )
        lay.addWidget(self.gap_label)

        self.dom_label = QLabel("Dominant: —")
        self.dom_label.setStyleSheet(f"color:{TEXT_MED}; font-size:11px;")
        self.dom_label.setWordWrap(True)
        self.tdm_label = QLabel("|TDM(S₁)| = —")
        self.tdm_label.setStyleSheet(f"color:{TEXT_MED}; font-size:11px;")
        self.p2_label = QLabel("2·P² = —")
        self.p2_label.setStyleSheet(f"color:{TEXT_MED}; font-size:11px;")
        lay.addWidget(self.dom_label)
        lay.addWidget(self.tdm_label)
        lay.addWidget(self.p2_label)

    def _make_state_box(self, title: str, tint: str) -> dict:
        frame = QFrame()
        frame.setStyleSheet(f"QFrame {{ background:{BG_CARD_HI}; "
                            f"border:1px solid {BORDER}; border-radius:4px; }}")
        v = QVBoxLayout(frame)
        v.setContentsMargins(8, 6, 8, 6)
        v.setSpacing(1)
        t = QLabel(title)
        t.setStyleSheet(f"color:{tint}; font-size:10px; font-weight:700; "
                        f"letter-spacing:1px;")
        energy = QLabel("— eV")
        energy.setStyleSheet(f"color:{TEXT_BRIGHT}; font-size:14px; "
                             f"font-weight:700;")
        wavelength = QLabel("—")
        wavelength.setStyleSheet(f"color:{TEXT_DIM}; font-size:10px;")
        v.addWidget(t)
        v.addWidget(energy)
        v.addWidget(wavelength)
        return {"frame": frame, "energy": energy, "wavelength": wavelength}

    def set_point(self, point: Optional[AnglePoint], threshold_ev: float):
        if point is None:
            self._clear()
            return
        s1, t1 = point.s1_state, point.t1_state
        if s1:
            self.s1_box["energy"].setText(f"{s1.energy_ev:.4f} eV")
            self.s1_box["wavelength"].setText(
                f"{s1.wavelength_nm:.1f} nm · f={s1.oscillator_strength:.3f}"
            )
        if t1:
            self.t1_box["energy"].setText(f"{t1.energy_ev:.4f} eV")
            self.t1_box["wavelength"].setText(
                f"{t1.wavelength_nm:.1f} nm · <S²>={t1.s_squared:.1f}"
            )

        gap = point.s1_t1_gap_ev
        if gap is not None:
            is_tadf = gap <= threshold_ev
            tint = GOOD if is_tadf else (WARN if gap <= 2 * threshold_ev else BAD)
            tag = " ✓ TADF" if is_tadf else ""
            self.gap_label.setText(f"ΔE = {gap:+.4f} eV{tag}")
            self.gap_label.setStyleSheet(
                f"color:{tint}; font-size:13px; font-weight:700; "
                f"padding:6px; border-radius:4px; background:{BG_CARD_HI};"
            )

        dom = point.s1_dominant_transition
        if dom:
            i, a, c = dom
            self.dom_label.setText(f"Dominant: MO {i}→{a} ({c:+.4f})")
        self.tdm_label.setText(f"|TDM(S₁)| = {point.s1_tdm_magnitude:.4f} Au")
        self.p2_label.setText(f"2·P² = {point.two_p_squared:.4f}")

    def _clear(self):
        self.s1_box["energy"].setText("— eV"); self.s1_box["wavelength"].setText("—")
        self.t1_box["energy"].setText("— eV"); self.t1_box["wavelength"].setText("—")
        self.gap_label.setText("ΔE = —")
        self.dom_label.setText("Dominant: —")
        self.tdm_label.setText("|TDM(S₁)| = —")
        self.p2_label.setText("2·P² = —")


# ── Main window ──────────────────────────────────────────────────────────────

class AngleScanWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Angle Scan — TADF Workbench")
        self.resize(1600, 1000)
        self._scan: Optional[AngleScan] = None
        self._current_angle: Optional[float] = None
        self._threshold = 0.2
        self.setStyleSheet(f"""
            QMainWindow, QWidget {{ background-color: {BG_WINDOW};
                                    color: {TEXT_BRIGHT}; }}
            QSplitter::handle {{ background: {BORDER}; }}
            QSplitter::handle:horizontal {{ width: 2px; }}
            QSplitter::handle:vertical   {{ height: 2px; }}
            QStatusBar {{ background: {BG_CARD}; color: {TEXT_MED};
                          border-top: 1px solid {BORDER}; }}
            QMenuBar {{ background: {BG_CARD}; color: {TEXT_BRIGHT};
                        border-bottom: 1px solid {BORDER}; padding: 3px; }}
            QMenuBar::item:selected {{ background: #094771; }}
            QMenu {{ background: {BG_CARD}; border: 1px solid {BORDER}; }}
            QMenu::item:selected {{ background: #094771; }}
        """)
        self._build_ui()
        self._build_menu()

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(8)

        # 1. Compact stat bar on top
        root.addWidget(self._build_stat_bar())

        # 2. Main splitter — MOLECULE (hero, left) | sidebar (right)
        main_split = QSplitter(Qt.Horizontal)

        # Left: 3D molecule viewer (dominant)
        self.viewer_3d = Dihedral3DPanel()
        self.viewer_3d.setMinimumWidth(640)
        viewer_frame = QFrame()
        viewer_frame.setStyleSheet(f"QFrame {{ background:{BG_CARD}; "
                                   f"border:1px solid {BORDER}; "
                                   f"border-radius:6px; }}")
        vf_lay = QVBoxLayout(viewer_frame)
        vf_lay.setContentsMargins(1, 1, 1, 1)
        vf_lay.setSpacing(0)
        vf_lay.addWidget(self.viewer_3d)
        main_split.addWidget(viewer_frame)

        # Right: narrow sidebar — current state + plots
        sidebar = QWidget()
        sidebar.setMinimumWidth(380)
        sidebar.setMaximumWidth(520)
        side_lay = QVBoxLayout(sidebar)
        side_lay.setContentsMargins(0, 0, 0, 0)
        side_lay.setSpacing(6)

        self.state_card = CurrentStateCard()
        self.state_card.setMaximumHeight(200)
        side_lay.addWidget(self.state_card)

        plots_split = QSplitter(Qt.Vertical)
        self.gap_plot = GapPlot()
        self.tdm_plot = TDMPlot()
        self.p2_plot = TwoPSquaredPlot()
        plots_split.addWidget(self.gap_plot)
        plots_split.addWidget(self.tdm_plot)
        plots_split.addWidget(self.p2_plot)
        plots_split.setStretchFactor(0, 3)  # gap plot biggest of the three
        plots_split.setStretchFactor(1, 2)
        plots_split.setStretchFactor(2, 2)
        side_lay.addWidget(plots_split, 1)

        main_split.addWidget(sidebar)
        main_split.setStretchFactor(0, 5)  # molecule dominant
        main_split.setStretchFactor(1, 2)
        root.addWidget(main_split, 1)

        # 3. Compact table at the bottom
        self.table = AngleScanTable()
        self.table.setMinimumHeight(150)
        self.table.setMaximumHeight(210)
        root.addWidget(self.table)

        # Cross-wiring
        self.table.angle_selected.connect(self._on_angle_selected)
        self.viewer_3d.angle_changed.connect(self._on_angle_selected)
        for p in (self.gap_plot, self.tdm_plot, self.p2_plot):
            p.angle_clicked.connect(self._on_angle_selected)

        self.statusBar().showMessage("No scan loaded.")

    def _build_stat_bar(self) -> QWidget:
        bar = QWidget()
        bar.setMaximumHeight(46)
        row = QHBoxLayout(bar)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        self.chip_theta = StatChip("CURRENT θ")
        self.chip_gap_now = StatChip("ΔE AT θ")
        self.chip_tadf = StatChip("TADF")
        self.chip_min = StatChip("MIN ΔE @ θ")
        self.chip_tdm = StatChip("|TDM|")
        self.chip_2p2 = StatChip("2·P²")
        for c in (self.chip_theta, self.chip_gap_now, self.chip_tadf,
                  self.chip_min, self.chip_tdm, self.chip_2p2):
            row.addWidget(c, 1)
        return bar

    def _build_menu(self):
        mb = self.menuBar()
        file_menu = mb.addMenu("&File")
        open_act = QAction("Open Scan Folder…", self,
                           shortcut=QKeySequence("Ctrl+Shift+O"))
        open_act.triggered.connect(self._prompt_open_folder)
        file_menu.addAction(open_act)
        import_act = QAction("Import Log Files…", self,
                             shortcut=QKeySequence("Ctrl+I"))
        import_act.triggered.connect(self._prompt_import_files)
        file_menu.addAction(import_act)
        file_menu.addSeparator()
        export_act = QAction("Export Summary CSV…", self)
        export_act.triggered.connect(self._export_csv)
        file_menu.addAction(export_act)

        scan_menu = mb.addMenu("&Scan")
        config_act = QAction("Configure Donor/Acceptor…", self)
        config_act.triggered.connect(self._configure_donor_acceptor)
        scan_menu.addAction(config_act)
        thr_act = QAction("Set ΔE Threshold (eV)…", self)
        thr_act.triggered.connect(self._set_threshold)
        scan_menu.addAction(thr_act)

    # ── Public API ──

    def set_scan(self, scan: AngleScan):
        self._scan = scan
        self._current_angle = None
        self.gap_plot.set_scan(scan)
        self.gap_plot.set_threshold(self._threshold)
        self.tdm_plot.set_scan(scan)
        self.p2_plot.set_scan(scan)
        self.table.set_scan(scan, threshold_ev=self._threshold)
        self.viewer_3d.set_scan(scan)
        if not scan.points:
            self.state_card.set_point(None, self._threshold)
            self._refresh_chips()
        self._update_status()

    # ── Handlers ──

    def _prompt_open_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Select Scan Folder")
        if not folder:
            return
        try:
            scan = parse_angle_folder(folder)
        except Exception as exc:
            QMessageBox.critical(self, "Parse Error", str(exc))
            return
        if not scan.points:
            QMessageBox.warning(self, "Empty Scan",
                                "No log files with detectable angles.")
            return
        self.set_scan(scan)
        if scan.donor_atom_index is None:
            self._configure_donor_acceptor()

    def _prompt_import_files(self):
        """Multi-select log files; merge into the current scan (or start one)."""
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Import Gaussian Log Files",
            "", "Gaussian logs (*.log *.out);;All Files (*)",
        )
        if not paths:
            return

        target = self._scan if self._scan is not None else AngleScan(name="Imported")
        added, replaced, skipped = 0, 0, 0
        skipped_details = []
        replace_all_remaining: Optional[bool] = None  # True/False once user chooses

        for path in paths:
            angle = angle_from_filename(path)
            if angle is None:
                # Prompt user for the angle
                value, ok = QInputDialog.getDouble(
                    self, "Assign Angle",
                    f"No angle in filename:\n{path}\n\nEnter θ (degrees):",
                    value=0.0, min=-360.0, max=360.0, decimals=2,
                )
                if not ok:
                    skipped += 1
                    skipped_details.append(f"{path} (no angle)")
                    continue
                angle = value

            try:
                molecule = parse_gaussian_file(path)
            except Exception as exc:  # noqa: BLE001
                skipped += 1
                skipped_details.append(f"{path} (parse error: {exc})")
                continue

            point = AnglePoint(angle_deg=angle, source_path=path, molecule=molecule)

            if target.find_point(angle) is not None:
                # Conflict — ask once, then remember for the rest of the batch
                if replace_all_remaining is None:
                    resp = QMessageBox.question(
                        self, "Angle Already in Scan",
                        f"θ = {angle:.2f}° is already present.\n\n"
                        "Replace it with the imported file?\n"
                        "(Yes = replace this and all later conflicts, "
                        "No = skip all later conflicts)",
                        QMessageBox.Yes | QMessageBox.No,
                    )
                    replace_all_remaining = (resp == QMessageBox.Yes)
                if not replace_all_remaining:
                    skipped += 1
                    skipped_details.append(f"{path} (θ {angle:.2f}° already present)")
                    continue

            was_replaced = target.add_or_replace(point)
            if was_replaced:
                replaced += 1
            else:
                added += 1

        if not added and not replaced:
            QMessageBox.information(
                self, "Nothing Imported",
                "No files imported.\n\n" + ("\n".join(skipped_details) if skipped_details else ""),
            )
            return

        # Adopt as current scan if we started fresh
        if self._scan is None:
            self.set_scan(target)
            if target.donor_atom_index is None:
                self._configure_donor_acceptor()
        else:
            # Refresh all views with the (mutated) existing scan
            self.set_scan(self._scan)

        # Status summary
        parts = []
        if added: parts.append(f"{added} added")
        if replaced: parts.append(f"{replaced} replaced")
        if skipped: parts.append(f"{skipped} skipped")
        self.statusBar().showMessage("Import: " + " · ".join(parts), 5000)

    def _configure_donor_acceptor(self):
        if not self._scan:
            QMessageBox.information(self, "No Scan", "Load a scan folder first.")
            return
        dlg = DonorAcceptorDialog(scan=self._scan, parent=self)
        if dlg.exec_():
            if self._current_angle is not None:
                self.viewer_3d.set_angle(self._current_angle)
            else:
                self.viewer_3d.set_scan(self._scan)

    def _set_threshold(self):
        val, ok = QInputDialog.getDouble(
            self, "ΔE Threshold", "TADF candidate threshold (eV):",
            value=self._threshold, min=0.0, max=2.0, decimals=3,
        )
        if ok:
            self._threshold = val
            if self._scan:
                self.gap_plot.set_threshold(val)
                self.table.set_scan(self._scan, threshold_ev=val)
                self._refresh_chips()
                self._update_status()

    def _export_csv(self):
        if not self._scan:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export Summary CSV", f"{self._scan.name}_summary.csv",
            "CSV (*.csv)")
        if not path:
            return
        import csv
        from ..analysis import summary_row
        keys = ["angle_deg", "s1_energy_ev", "t1_energy_ev", "gap_ev",
                "tdm_magnitude", "two_p_squared", "s1_dominant", "source_path"]
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(keys)
            for p in self._scan.sorted_points:
                row = summary_row(p)
                w.writerow([row.get(k) for k in keys])
        self.statusBar().showMessage(f"Exported: {path}", 4000)

    def _on_angle_selected(self, angle_deg: float):
        if self._current_angle == angle_deg:
            return
        self._current_angle = angle_deg
        self.viewer_3d.set_angle(angle_deg)
        if self._scan:
            pts = self._scan.sorted_points
            for row, p in enumerate(pts):
                if p.angle_deg == angle_deg:
                    self.table.blockSignals(True)
                    self.table.selectRow(row)
                    self.table.blockSignals(False)
                    break
        self._refresh_chips()
        self._update_status()

    # ── Rendering helpers ──

    def _refresh_chips(self):
        scan = self._scan
        if not scan or not scan.points:
            for chip in (self.chip_theta, self.chip_gap_now, self.chip_tadf,
                         self.chip_min, self.chip_tdm, self.chip_2p2):
                chip.set_value("—")
            self.state_card.set_point(None, self._threshold)
            return

        point = None
        if self._current_angle is not None:
            point = next(
                (p for p in scan.points if p.angle_deg == self._current_angle),
                None,
            )
        self.state_card.set_point(point, self._threshold)

        if point:
            self.chip_theta.set_value(f"{point.angle_deg:.0f}°")
            gap = point.s1_t1_gap_ev
            if gap is not None:
                tint = (GOOD if gap <= self._threshold
                        else WARN if gap <= 2 * self._threshold else BAD)
                self.chip_gap_now.set_value(f"{gap:+.3f} eV", "S₁−T₁", tint=tint)
                is_tadf = gap <= self._threshold
                self.chip_tadf.set_value(
                    "YES" if is_tadf else "no",
                    f"≤ {self._threshold:.2f} eV",
                    tint=GOOD if is_tadf else TEXT_DIM,
                )
            else:
                self.chip_gap_now.set_value("—")
                self.chip_tadf.set_value("—")
            self.chip_tdm.set_value(f"{point.s1_tdm_magnitude:.3f}", "Au")
            self.chip_2p2.set_value(f"{point.two_p_squared:.3f}", "S₁ CI")
        else:
            for chip in (self.chip_theta, self.chip_gap_now, self.chip_tadf,
                         self.chip_tdm, self.chip_2p2):
                chip.set_value("—")

        best = scan.minimum_gap_point
        if best and best.s1_t1_gap_ev is not None:
            tint = GOOD if best.s1_t1_gap_ev <= self._threshold else WARN
            self.chip_min.set_value(
                f"{best.s1_t1_gap_ev:+.3f} eV",
                f"@ {best.angle_deg:.0f}°", tint=tint,
            )
        else:
            self.chip_min.set_value("—")

    def _update_status(self):
        if not self._scan:
            self.statusBar().showMessage("No scan loaded.")
            return
        msg = f"{len(self._scan)} pts"
        best = self._scan.minimum_gap_point
        if best is not None and best.s1_t1_gap_ev is not None:
            msg += f" | min gap {best.s1_t1_gap_ev:.4f} eV @ {best.angle_deg:.1f}°"
        msg += f" | threshold {self._threshold:.3f} eV"
        n_tadf = len(self._scan.points_below_threshold(self._threshold))
        msg += f" | {n_tadf} TADF candidates"
        self.statusBar().showMessage(msg)
