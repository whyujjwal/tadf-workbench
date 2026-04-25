"""
Basic tests for molecule visualizer core components.
"""

import pytest
import numpy as np
from tadf_workbench.core import Atom, Stage, Molecule
from tadf_workbench.utils import get_element_symbol, get_atom_color


class TestAtom:
    """Test Atom class."""
    
    def test_atom_creation(self):
        atom = Atom(6, 1.0, 2.0, 3.0, 0)
        assert atom.atomic_number == 6
        assert atom.x == 1.0
        assert atom.y == 2.0
        assert atom.z == 3.0
    
    def test_atom_position(self):
        atom = Atom(6, 1.0, 2.0, 3.0, 0)
        pos = atom.position
        assert isinstance(pos, np.ndarray)
        assert np.allclose(pos, [1.0, 2.0, 3.0])
    
    def test_distance_between_atoms(self):
        atom1 = Atom(6, 0.0, 0.0, 0.0, 0)
        atom2 = Atom(6, 3.0, 4.0, 0.0, 1)
        distance = atom1.distance_to(atom2)
        assert np.isclose(distance, 5.0)


class TestStage:
    """Test Stage class."""
    
    def test_stage_creation(self):
        stage = Stage("Test Stage")
        assert stage.name == "Test Stage"
        assert stage.atom_count == 0
    
    def test_add_atom(self):
        stage = Stage("Test")
        stage.add_atom(6, 0.0, 0.0, 0.0)
        assert stage.atom_count == 1
        assert stage.atoms[0].atomic_number == 6
    
    def test_stage_from_data(self):
        atoms = [(6, 0.0, 0.0, 0.0), (1, 1.0, 0.0, 0.0)]
        stage = Stage("Test", atoms)
        assert stage.atom_count == 2
    
    def test_center_of_mass(self):
        atoms = [
            (6, 0.0, 0.0, 0.0),
            (6, 2.0, 0.0, 0.0),
            (6, 0.0, 2.0, 0.0),
        ]
        stage = Stage("Test", atoms)
        center = stage.get_center_of_mass()
        expected = np.array([2.0/3, 2.0/3, 0.0])
        assert np.allclose(center, expected)
    
    def test_find_bonds(self):
        # Two carbon atoms close together
        atoms = [
            (6, 0.0, 0.0, 0.0),
            (6, 1.5, 0.0, 0.0),  # Within bond distance
            (6, 5.0, 0.0, 0.0),  # Too far
        ]
        stage = Stage("Test", atoms)
        bonds = stage.find_bonds()
        assert len(bonds) == 1
        assert bonds[0] == (0, 1)


class TestMolecule:
    """Test Molecule class."""
    
    def test_molecule_creation(self):
        mol = Molecule("Test Molecule")
        assert mol.name == "Test Molecule"
        assert mol.stage_count == 0
    
    def test_add_stage(self):
        mol = Molecule("Test")
        stage = Stage("Stage 1")
        mol.add_stage(stage)
        assert mol.stage_count == 1
    
    def test_navigation(self):
        mol = Molecule("Test")
        mol.add_stage(Stage("Stage 1"))
        mol.add_stage(Stage("Stage 2"))
        mol.add_stage(Stage("Stage 3"))
        
        assert mol.current_stage_index == 0
        assert mol.next_stage() is True
        assert mol.current_stage_index == 1
        assert mol.previous_stage() is True
        assert mol.current_stage_index == 0
        assert mol.previous_stage() is False  # Can't go before first
    
    def test_to_dict(self):
        mol = Molecule("Test")
        atoms = [(6, 0.0, 0.0, 0.0)]
        mol.add_stage_from_data("Stage 1", atoms)
        
        data = mol.to_dict()
        assert data['name'] == "Test"
        assert len(data['stages']) == 1
        assert data['stages'][0]['name'] == "Stage 1"


class TestUtils:
    """Test utility functions."""
    
    def test_get_element_symbol(self):
        assert get_element_symbol(1) == 'H'
        assert get_element_symbol(6) == 'C'
        assert get_element_symbol(8) == 'O'
        assert get_element_symbol(999) == '#999'
    
    def test_get_atom_color(self):
        color = get_atom_color(1)  # Hydrogen
        assert len(color) == 4  # RGBA
        assert color == (1.0, 1.0, 1.0, 1.0)  # White
        
        color = get_atom_color(8)  # Oxygen
        assert color == (1.0, 0.0, 0.0, 1.0)  # Red


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
