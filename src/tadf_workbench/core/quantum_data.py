"""
Quantum chemistry data structures.
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional
import numpy as np


@dataclass
class DipoleMoment:
    """Electric dipole moment."""
    x: float
    y: float
    z: float
    total: float
    
    @property
    def vector(self) -> np.ndarray:
        """Get dipole as vector."""
        return np.array([self.x, self.y, self.z])
    
    @property
    def magnitude(self) -> float:
        """Get dipole magnitude."""
        return self.total
    
    def __str__(self):
        return f"μ = {self.total:.4f} D [{self.x:.4f}, {self.y:.4f}, {self.z:.4f}]"


@dataclass
class QuadrupoleMoment:
    """Electric quadrupole moment tensor."""
    xx: float
    yy: float
    zz: float
    xy: float
    xz: float
    yz: float
    
    @property
    def tensor(self) -> np.ndarray:
        """Get quadrupole as symmetric tensor."""
        return np.array([
            [self.xx, self.xy, self.xz],
            [self.xy, self.yy, self.yz],
            [self.xz, self.yz, self.zz]
        ])


@dataclass
class ExcitedState:
    """Excited state information."""
    number: int
    multiplicity: str  # 'Singlet' or 'Triplet'
    symmetry: str
    energy_ev: float
    wavelength_nm: float
    oscillator_strength: float
    s_squared: float
    
    # Transition vectors
    electric_dipole: Optional[Tuple[float, float, float]] = None
    magnetic_dipole: Optional[Tuple[float, float, float]] = None
    velocity_dipole: Optional[Tuple[float, float, float]] = None
    
    # Additional properties
    rotatory_strength: Optional[float] = None
    orbital_transitions: List[Tuple[int, int, float]] = None
    
    @property
    def is_singlet(self) -> bool:
        """Check if state is singlet."""
        return 'Singlet' in self.multiplicity
    
    @property
    def is_triplet(self) -> bool:
        """Check if state is triplet."""
        return 'Triplet' in self.multiplicity
    
    @property
    def transition_strength(self) -> float:
        """Get transition dipole magnitude."""
        if self.electric_dipole:
            x, y, z = self.electric_dipole
            return np.sqrt(x**2 + y**2 + z**2)
        return 0.0
    
    def __str__(self):
        return (f"State {self.number}: {self.multiplicity}-{self.symmetry} "
                f"{self.energy_ev:.4f} eV ({self.wavelength_nm:.2f} nm) f={self.oscillator_strength:.4f}")


@dataclass
class TransitionDipole:
    """Ground to excited state transition electric dipole moment."""
    state: int
    x: float  # X component (Au)
    y: float  # Y component (Au)
    z: float  # Z component (Au)
    dip_strength: float  # Dipole strength (Dip. S.)
    osc_strength: float  # Oscillator strength (Osc.)

    # Velocity gauge dipole
    vel_x: float = 0.0
    vel_y: float = 0.0
    vel_z: float = 0.0
    vel_dip_strength: float = 0.0
    vel_osc_strength: float = 0.0

    # Magnetic dipole
    mag_x: float = 0.0
    mag_y: float = 0.0
    mag_z: float = 0.0

    # Rotatory strengths (cgs units)
    rotatory_velocity: float = 0.0
    rotatory_length: float = 0.0

    @property
    def vector(self) -> np.ndarray:
        """Get transition dipole as vector."""
        return np.array([self.x, self.y, self.z])

    @property
    def magnitude(self) -> float:
        """Get transition dipole magnitude."""
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def magnetic_magnitude(self) -> float:
        """Get magnetic dipole magnitude."""
        return np.sqrt(self.mag_x**2 + self.mag_y**2 + self.mag_z**2)

    @property
    def dominant_axis(self) -> str:
        """Get the dominant axis of the transition."""
        abs_vals = [abs(self.x), abs(self.y), abs(self.z)]
        max_idx = abs_vals.index(max(abs_vals))
        return ['X', 'Y', 'Z'][max_idx]

    def __str__(self):
        return f"State {self.state}: X={self.x:.4f}, Y={self.y:.4f}, Z={self.z:.4f}, Dip.S.={self.dip_strength:.4f}, Osc.={self.osc_strength:.4f}"


@dataclass
class MullikenCharges:
    """Mulliken atomic charges."""
    atom_indices: List[int]
    charges: List[float]
    
    def get_charge(self, atom_index: int) -> float:
        """Get charge for specific atom."""
        if atom_index < len(self.charges):
            return self.charges[atom_index]
        return 0.0
    
    @property
    def min_charge(self) -> float:
        """Get minimum charge."""
        return min(self.charges) if self.charges else 0.0
    
    @property
    def max_charge(self) -> float:
        """Get maximum charge."""
        return max(self.charges) if self.charges else 0.0


@dataclass
class MolecularOrbitals:
    """Molecular orbital energies."""
    energies: List[float]
    occupations: List[int]
    
    @property
    def homo_index(self) -> int:
        """Get HOMO index."""
        for i in range(len(self.occupations) - 1, -1, -1):
            if self.occupations[i] > 0:
                return i
        return -1
    
    @property
    def lumo_index(self) -> int:
        """Get LUMO index."""
        homo = self.homo_index
        if homo >= 0 and homo < len(self.energies) - 1:
            return homo + 1
        return -1
    
    @property
    def homo_energy(self) -> Optional[float]:
        """Get HOMO energy."""
        idx = self.homo_index
        return self.energies[idx] if idx >= 0 else None
    
    @property
    def lumo_energy(self) -> Optional[float]:
        """Get LUMO energy."""
        idx = self.lumo_index
        return self.energies[idx] if idx >= 0 else None
    
    @property
    def homo_lumo_gap(self) -> Optional[float]:
        """Get HOMO-LUMO gap in eV."""
        homo = self.homo_energy
        lumo = self.lumo_energy
        if homo is not None and lumo is not None:
            return lumo - homo
        return None


class QuantumData:
    """Container for all quantum chemistry data."""

    def __init__(self):
        self.dipole_moment: Optional[DipoleMoment] = None
        self.quadrupole_moment: Optional[QuadrupoleMoment] = None
        self.excited_states: List[ExcitedState] = []
        self.transition_dipoles: List[TransitionDipole] = []  # Transition electric dipole moments
        self.mulliken_charges: Optional[MullikenCharges] = None
        self.molecular_orbitals: Optional[MolecularOrbitals] = None

        # Octapole and hexadecapole (for advanced users)
        self.octapole_moment: Optional[dict] = None
        self.hexadecapole_moment: Optional[dict] = None

    def add_excited_state(self, state: ExcitedState):
        """Add an excited state."""
        self.excited_states.append(state)

    def add_transition_dipole(self, td: TransitionDipole):
        """Add a transition dipole moment."""
        self.transition_dipoles.append(td)

    def get_transition_dipoles_sorted(self, by: str = 'state', reverse: bool = False) -> List[TransitionDipole]:
        """Get transition dipoles sorted by specified field.

        Args:
            by: Sort field - 'state', 'dip_strength', 'osc_strength', 'x', 'y', 'z'
            reverse: If True, sort in descending order
        """
        key_map = {
            'state': lambda td: td.state,
            'dip_strength': lambda td: td.dip_strength,
            'osc_strength': lambda td: td.osc_strength,
            'x': lambda td: abs(td.x),
            'y': lambda td: abs(td.y),
            'z': lambda td: abs(td.z),
        }
        key_func = key_map.get(by, key_map['state'])
        return sorted(self.transition_dipoles, key=key_func, reverse=reverse)
    
    @property
    def num_excited_states(self) -> int:
        """Get number of excited states."""
        return len(self.excited_states)
    
    @property
    def singlet_states(self) -> List[ExcitedState]:
        """Get all singlet states."""
        return [s for s in self.excited_states if s.is_singlet]
    
    @property
    def triplet_states(self) -> List[ExcitedState]:
        """Get all triplet states."""
        return [s for s in self.excited_states if s.is_triplet]
    
    @property
    def strongest_transition(self) -> Optional[ExcitedState]:
        """Get state with strongest oscillator strength."""
        if not self.excited_states:
            return None
        return max(self.excited_states, key=lambda s: s.oscillator_strength)
    
    def get_state(self, number: int) -> Optional[ExcitedState]:
        """Get excited state by number."""
        for state in self.excited_states:
            if state.number == number:
                return state
        return None
    
    def to_dict(self) -> dict:
        """Convert to dictionary for serialization."""
        return {
            'dipole_moment': {
                'x': self.dipole_moment.x,
                'y': self.dipole_moment.y,
                'z': self.dipole_moment.z,
                'total': self.dipole_moment.total
            } if self.dipole_moment else None,
            'num_excited_states': self.num_excited_states,
            'excited_states': [
                {
                    'number': s.number,
                    'multiplicity': s.multiplicity,
                    'energy_ev': s.energy_ev,
                    'wavelength_nm': s.wavelength_nm,
                    'oscillator_strength': s.oscillator_strength
                }
                for s in self.excited_states
            ]
        }
