"""Batch-parse a folder of Gaussian log files representing a dihedral angle scan."""

import os
import re
from typing import Optional, List

from ..core import AngleScan, AnglePoint
from .gaussian_parser import parse_gaussian_file


_ANGLE_PATTERN = re.compile(r"(\d+)")

_SUPPORTED_EXTS = {".log", ".out"}


def angle_from_filename(path: str) -> Optional[float]:
    """Extract the dihedral angle (degrees) from the filename stem.

    Strategy: take the basename (no directory), drop extension, return the first
    integer substring. Returns None if no integer is present.
    """
    basename = os.path.basename(path)
    stem, _ = os.path.splitext(basename)
    match = _ANGLE_PATTERN.search(stem)
    if not match:
        return None
    return float(match.group(1))


def parse_angle_folder(folder: str) -> AngleScan:
    """Parse every supported log in `folder`, attaching its dihedral angle.

    Files without a detectable angle in the filename are skipped with a warning.
    Files that fail to parse are skipped with a warning.
    The returned AngleScan is named after the folder basename.
    """
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"Not a directory: {folder}")

    points: List[AnglePoint] = []
    for entry in sorted(os.listdir(folder)):
        full = os.path.join(folder, entry)
        if not os.path.isfile(full):
            continue
        _, ext = os.path.splitext(entry)
        if ext.lower() not in _SUPPORTED_EXTS:
            continue
        angle = angle_from_filename(entry)
        if angle is None:
            print(f"[angle_scan] Skipping '{entry}': no angle in filename")
            continue
        try:
            molecule = parse_gaussian_file(full)
        except Exception as exc:  # noqa: BLE001 — keep scan going
            print(f"[angle_scan] Skipping '{entry}': parse error: {exc}")
            continue
        points.append(AnglePoint(
            angle_deg=angle, source_path=full, molecule=molecule,
        ))

    name = os.path.basename(os.path.normpath(folder)) or "angle_scan"
    return AngleScan(name=name, points=points)
