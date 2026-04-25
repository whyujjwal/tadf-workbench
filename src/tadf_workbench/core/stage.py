"""
Core data model for representing a molecular stage.
"""

from typing import List, Tuple
from dataclasses import dataclass
import numpy as np


@dataclass
class Atom:
    """Represents a single atom with its properties."""
    
    atomic_number: int
    x: float
    y: float
    z: float
    index: int = 0
    
    @property
    def position(self) -> np.ndarray:
        """Get atom position as numpy array."""
        return np.array([self.x, self.y, self.z])
    
    def distance_to(self, other: 'Atom') -> float:
        """Calculate distance to another atom."""
        return np.linalg.norm(self.position - other.position)
    
    def __repr__(self) -> str:
        from ..utils import get_element_symbol
        symbol = get_element_symbol(self.atomic_number)
        return f"Atom({symbol}, x={self.x:.3f}, y={self.y:.3f}, z={self.z:.3f})"


class Stage:
    """Represents a single stage of molecular coordinates."""
    
    def __init__(self, name: str, atoms: List[Tuple[int, float, float, float]] = None):
        """
        Initialize a stage.
        
        Args:
            name: Name/description of this stage
            atoms: List of (atomic_number, x, y, z) tuples
        """
        self.name = name
        self._atoms: List[Atom] = []
        
        if atoms:
            for idx, (atomic_num, x, y, z) in enumerate(atoms):
                self.add_atom(atomic_num, x, y, z, idx)
    
    def add_atom(self, atomic_number: int, x: float, y: float, z: float, index: int = None):
        """Add an atom to this stage."""
        if index is None:
            index = len(self._atoms)
        atom = Atom(atomic_number, x, y, z, index)
        self._atoms.append(atom)
    
    @property
    def atoms(self) -> List[Atom]:
        """Get list of atoms in this stage."""
        return self._atoms
    
    @property
    def atom_count(self) -> int:
        """Get number of atoms in this stage."""
        return len(self._atoms)
    
    @property
    def positions(self) -> np.ndarray:
        """Get all atom positions as numpy array (N x 3)."""
        return np.array([atom.position for atom in self._atoms])
    
    @property
    def atomic_numbers(self) -> List[int]:
        """Get list of atomic numbers."""
        return [atom.atomic_number for atom in self._atoms]
    
    def get_center_of_mass(self) -> np.ndarray:
        """Calculate geometric center (centroid) of all atoms."""
        if not self._atoms:
            return np.array([0.0, 0.0, 0.0])
        return np.mean(self.positions, axis=0)
    
    def get_bounds(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get bounding box of the molecule.
        
        Returns:
            Tuple of (min_coords, max_coords) as numpy arrays
        """
        if not self._atoms:
            return np.array([0, 0, 0]), np.array([0, 0, 0])
        
        positions = self.positions
        return np.min(positions, axis=0), np.max(positions, axis=0)
    
    def find_bonds(self, tolerance: float = 0.3) -> List[Tuple[int, int]]:
        """
        Find bonds between atoms based on distance.
        
        Args:
            tolerance: Additional tolerance for bond detection (Angstroms)
            
        Returns:
            List of (atom1_index, atom2_index) tuples representing bonds
        """
        from ..utils import calculate_bond_threshold
        
        bonds = []
        for i in range(len(self._atoms)):
            for j in range(i + 1, len(self._atoms)):
                atom1 = self._atoms[i]
                atom2 = self._atoms[j]
                
                threshold = calculate_bond_threshold(
                    atom1.atomic_number,
                    atom2.atomic_number,
                    tolerance
                )
                
                distance = atom1.distance_to(atom2)
                if distance < threshold:
                    bonds.append((i, j))
        
        return bonds
    
    def translate(self, offset: np.ndarray):
        """
        Translate all atoms by the given offset.
        
        Args:
            offset: 3D offset vector [dx, dy, dz]
        """
        for atom in self._atoms:
            atom.x += offset[0]
            atom.y += offset[1]
            atom.z += offset[2]
    
    def center_at_origin(self):
        """Center the molecule at the origin."""
        center = self.get_center_of_mass()
        self.translate(-center)
    
    def to_dict(self) -> dict:
        """Convert stage to dictionary format."""
        return {
            'name': self.name,
            'atoms': [
                (atom.atomic_number, atom.x, atom.y, atom.z)
                for atom in self._atoms
            ]
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> 'Stage':
        """Create stage from dictionary format."""
        return cls(name=data['name'], atoms=data.get('atoms', []))
    
    def __repr__(self) -> str:
        return f"Stage(name='{self.name}', atoms={self.atom_count})"
    
    def __len__(self) -> int:
        return self.atom_count
