"""
Molecule data model representing multiple stages.
"""

from typing import List, Optional
from .stage import Stage


class Molecule:
    """Represents a molecule with multiple stages/conformations."""
    
    def __init__(self, name: str = "Unnamed"):
        """
        Initialize a new molecule.
        
        Args:
            name: Name of the molecule
        """
        self.name = name
        self._stages: List[Stage] = []
        self._current_stage_index = 0
        self.quantum_data = None  # Will hold QuantumData object
    
    def add_stage(self, stage: Stage):
        """Add a stage to the molecule."""
        self._stages.append(stage)
    
    def add_stage_from_data(self, name: str, atoms: list):
        """
        Add a stage from raw data.
        
        Args:
            name: Stage name
            atoms: List of (atomic_number, x, y, z) tuples
        """
        stage = Stage(name, atoms)
        self.add_stage(stage)
    
    @property
    def stages(self) -> List[Stage]:
        """Get all stages."""
        return self._stages
    
    @property
    def stage_count(self) -> int:
        """Get number of stages."""
        return len(self._stages)
    
    @property
    def current_stage(self) -> Optional[Stage]:
        """Get current stage."""
        if not self._stages:
            return None
        return self._stages[self._current_stage_index]
    
    @property
    def current_stage_index(self) -> int:
        """Get current stage index."""
        return self._current_stage_index
    
    @current_stage_index.setter
    def current_stage_index(self, index: int):
        """Set current stage index."""
        if 0 <= index < len(self._stages):
            self._current_stage_index = index
    
    def get_stage(self, index: int) -> Optional[Stage]:
        """Get stage by index."""
        if 0 <= index < len(self._stages):
            return self._stages[index]
        return None
    
    def next_stage(self) -> bool:
        """
        Move to next stage.
        
        Returns:
            True if moved, False if already at last stage
        """
        if self._current_stage_index < len(self._stages) - 1:
            self._current_stage_index += 1
            return True
        return False
    
    def previous_stage(self) -> bool:
        """
        Move to previous stage.
        
        Returns:
            True if moved, False if already at first stage
        """
        if self._current_stage_index > 0:
            self._current_stage_index -= 1
            return True
        return False
    
    def clear_stages(self):
        """Remove all stages."""
        self._stages.clear()
        self._current_stage_index = 0
    
    def get_stage_names(self) -> List[str]:
        """Get list of all stage names."""
        return [stage.name for stage in self._stages]
    
    def to_dict(self) -> dict:
        """Convert molecule to dictionary format."""
        return {
            'name': self.name,
            'current_stage': self._current_stage_index,
            'stages': [stage.to_dict() for stage in self._stages]
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> 'Molecule':
        """Create molecule from dictionary format."""
        molecule = cls(name=data.get('name', 'Untitled Molecule'))
        for stage_data in data.get('stages', []):
            stage = Stage.from_dict(stage_data)
            molecule.add_stage(stage)
        molecule.current_stage_index = data.get('current_stage', 0)
        return molecule
    
    def __repr__(self) -> str:
        return f"Molecule(name='{self.name}', stages={self.stage_count})"
    
    def __len__(self) -> int:
        return self.stage_count
