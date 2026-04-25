"""Entry point: launch the TADF angle-scan workbench."""
import sys

from PyQt5.QtWidgets import QApplication

from .ui import AngleScanWindow


def main() -> int:
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setApplicationName("TADF Workbench")
    app.setApplicationVersion("0.1.0")
    win = AngleScanWindow()
    win.show()
    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())
