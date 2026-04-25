"""
Atom properties and utilities for molecular visualization.
"""

from typing import Tuple


# Atomic numbers to element symbols
ATOMIC_NUMBERS = {
    1: 'H',   # Hydrogen
    6: 'C',   # Carbon
    7: 'N',   # Nitrogen
    8: 'O',   # Oxygen
    9: 'F',   # Fluorine
    15: 'P',  # Phosphorus
    16: 'S',  # Sulfur
    17: 'Cl', # Chlorine
    35: 'Br', # Bromine
    53: 'I',  # Iodine
}

# Element symbols to atomic numbers
ELEMENT_SYMBOLS = {v: k for k, v in ATOMIC_NUMBERS.items()}


def get_element_symbol(atomic_number: int) -> str:
    """
    Get element symbol from atomic number.
    
    Args:
        atomic_number: Atomic number of the element
        
    Returns:
        Element symbol (e.g., 'C', 'H', 'O')
    """
    return ATOMIC_NUMBERS.get(atomic_number, f'#{atomic_number}')


def get_atomic_number(symbol: str) -> int:
    """
    Get atomic number from element symbol.
    
    Args:
        symbol: Element symbol (e.g., 'C', 'H', 'O')
        
    Returns:
        Atomic number
    """
    return ELEMENT_SYMBOLS.get(symbol.upper(), 0)


def get_atom_color(atomic_number: int) -> Tuple[float, float, float, float]:
    """
    Get RGBA color for atom visualization based on CPK coloring scheme.
    
    Args:
        atomic_number: Atomic number of the element
        
    Returns:
        RGBA tuple (r, g, b, a) with values from 0.0 to 1.0
    """
    colors = {
        1: (1.0, 1.0, 1.0, 1.0),    # H - White
        6: (0.3, 0.3, 0.3, 1.0),    # C - Dark Gray
        7: (0.2, 0.4, 1.0, 1.0),    # N - Blue
        8: (1.0, 0.0, 0.0, 1.0),    # O - Red
        9: (0.0, 1.0, 0.5, 1.0),    # F - Light Green
        15: (1.0, 0.5, 0.0, 1.0),   # P - Orange
        16: (1.0, 1.0, 0.0, 1.0),   # S - Yellow
        17: (0.0, 1.0, 0.0, 1.0),   # Cl - Green
        35: (0.6, 0.1, 0.1, 1.0),   # Br - Dark Red
        53: (0.5, 0.0, 0.5, 1.0),   # I - Purple
    }
    return colors.get(atomic_number, (0.5, 0.5, 0.5, 1.0))  # Default: Gray


def get_atom_radius(atomic_number: int) -> float:
    """
    Get Van der Waals radius for atom visualization (in Angstroms).
    
    Args:
        atomic_number: Atomic number of the element
        
    Returns:
        Radius in Angstroms (scaled for visualization)
    """
    radii = {
        1: 0.3,   # H
        6: 0.5,   # C
        7: 0.5,   # N
        8: 0.5,   # O
        9: 0.4,   # F
        15: 0.6,  # P
        16: 0.6,  # S
        17: 0.55, # Cl
        35: 0.65, # Br
        53: 0.75, # I
    }
    return radii.get(atomic_number, 0.4)


def get_covalent_radius(atomic_number: int) -> float:
    """
    Get covalent radius for bond distance calculations (in Angstroms).
    
    Args:
        atomic_number: Atomic number of the element
        
    Returns:
        Covalent radius in Angstroms
    """
    radii = {
        1: 0.31,   # H
        6: 0.76,   # C
        7: 0.71,   # N
        8: 0.66,   # O
        9: 0.57,   # F
        15: 1.07,  # P
        16: 1.05,  # S
        17: 1.02,  # Cl
        35: 1.20,  # Br
        53: 1.39,  # I
    }
    return radii.get(atomic_number, 0.7)


def calculate_bond_threshold(atom1_num: int, atom2_num: int, tolerance: float = 0.3) -> float:
    """
    Calculate maximum distance for bond detection between two atoms.
    
    Args:
        atom1_num: Atomic number of first atom
        atom2_num: Atomic number of second atom
        tolerance: Additional tolerance factor (default: 0.3 Angstroms)
        
    Returns:
        Maximum bond distance in Angstroms
    """
    r1 = get_covalent_radius(atom1_num)
    r2 = get_covalent_radius(atom2_num)
    return r1 + r2 + tolerance


def get_element_name(atomic_number: int) -> str:
    """
    Get full element name from atomic number.
    
    Args:
        atomic_number: Atomic number of the element
        
    Returns:
        Full element name
    """
    names = {
        1: 'Hydrogen',
        6: 'Carbon',
        7: 'Nitrogen',
        8: 'Oxygen',
        9: 'Fluorine',
        15: 'Phosphorus',
        16: 'Sulfur',
        17: 'Chlorine',
        35: 'Bromine',
        53: 'Iodine',
    }
    return names.get(atomic_number, f'Element {atomic_number}')
