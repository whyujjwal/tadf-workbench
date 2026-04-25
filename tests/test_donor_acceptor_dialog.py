# tests/test_donor_acceptor_dialog.py
import pytest
pytest.importorskip("PyQt5")
from PyQt5.QtWidgets import QDialog
from tadf_workbench.ui.donor_acceptor_dialog import DonorAcceptorDialog
from tadf_workbench.parsers import parse_angle_folder
from tests.fixtures.build_fixtures import write_scan_folder


@pytest.fixture
def scan(tmp_path):
    return parse_angle_folder(str(write_scan_folder(tmp_path / "s", {10: (3.4, 2.9)})))


def test_dialog_requires_selections_before_accept(qtbot, scan):
    dlg = DonorAcceptorDialog(scan=scan)
    qtbot.addWidget(dlg)
    assert not dlg.ok_button.isEnabled()
    dlg.set_selection(donor=0, acceptor=1, bond=(0, 1))
    assert dlg.ok_button.isEnabled()


def test_dialog_rejects_bond_with_equal_atoms(qtbot, scan):
    dlg = DonorAcceptorDialog(scan=scan)
    qtbot.addWidget(dlg)
    dlg.set_selection(donor=0, acceptor=1, bond=(0, 0))
    assert not dlg.ok_button.isEnabled()


def test_dialog_writes_back_to_scan_on_accept(qtbot, scan):
    dlg = DonorAcceptorDialog(scan=scan)
    qtbot.addWidget(dlg)
    dlg.set_selection(donor=0, acceptor=1, bond=(0, 1))
    dlg.accept()
    assert scan.donor_atom_index == 0
    assert scan.acceptor_atom_index == 1
    assert scan.rotatable_bond == (0, 1)
