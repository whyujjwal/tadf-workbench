"""Convenience launcher: `python run.py [scan_folder]`.

If a folder path is supplied, it is opened immediately as an angle scan and
the donor/acceptor configuration dialog is shown.
"""
import sys
from pathlib import Path

from PyQt5.QtWidgets import QApplication

from tadf_workbench.parsers import parse_angle_folder
from tadf_workbench.ui import AngleScanWindow, DonorAcceptorDialog


def main() -> int:
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    win = AngleScanWindow()
    if len(sys.argv) > 1:
        folder = Path(sys.argv[1]).expanduser().resolve()
        if not folder.is_dir():
            print(f"not a directory: {folder}", file=sys.stderr)
            return 2
        scan = parse_angle_folder(str(folder))
        if scan.points:
            win.set_scan(scan)
            DonorAcceptorDialog(scan=scan, parent=win).exec_()
            win.set_scan(scan)  # refresh now that fragments are configured
    win.show()
    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())
