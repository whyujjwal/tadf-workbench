"""Convenience launcher: `python run.py [scan_folder]`.

If a folder path is supplied, it is opened immediately as an angle scan and
the donor/acceptor configuration dialog is shown.
"""
import sys
from pathlib import Path

# Make the package importable from a fresh checkout (no editable install needed)
_SRC = Path(__file__).resolve().parent / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

try:
    from PyQt5.QtWidgets import QApplication
except ModuleNotFoundError:
    sys.stderr.write(
        "PyQt5 is not installed in this Python environment.\n"
        "Activate the project venv first, e.g.:\n\n"
        "    source .venv/bin/activate\n"
        "    python run.py data/demo_scan\n\n"
        "Or run with the venv's Python directly:\n\n"
        "    .venv/bin/python run.py data/demo_scan\n"
    )
    sys.exit(1)

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
