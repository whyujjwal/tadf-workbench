# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller build spec for the TADF Workbench (one-folder build).

Produces a self-contained folder with the launcher executable plus all the
Qt/OpenGL/numpy runtime it needs — no Python install required on the target
machine.

    pyinstaller tadf-workbench.spec --noconfirm

Output (Windows):  dist/TADF-Workbench/TADF-Workbench.exe  (+ supporting files)

This same spec is used by the GitHub Actions Windows build
(.github/workflows/build-windows.yml). It also works locally on Windows.
"""
from PyInstaller.utils.hooks import collect_submodules, collect_data_files

# PyOpenGL loads its platform backend and array-format handlers dynamically via
# ctypes, so PyInstaller's static import analysis misses them and the 3D viewer
# crashes at runtime ("Unable to load array format handler"). Collecting every
# submodule is the documented, reliable fix. We do the same for pyqtgraph, which
# imports many of its drawing/template modules lazily.
hiddenimports = collect_submodules("OpenGL") + collect_submodules("pyqtgraph")

# pyqtgraph ships icons / CSS / .ui template data alongside its code.
datas = collect_data_files("pyqtgraph")

a = Analysis(
    ["run.py"],
    pathex=["src"],          # so `tadf_workbench` resolves during analysis
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        # Pin pyqtgraph to the PyQt5 binding and keep the bundle lean: if other
        # Qt bindings get pulled in, pyqtgraph may pick the wrong one at runtime.
        "PyQt6", "PySide2", "PySide6",
        "PyQt5.QtWebEngine", "PyQt5.QtWebEngineWidgets", "PyQt5.QtBluetooth",
        # Heavy libraries this app never imports.
        "tkinter", "matplotlib", "pytest", "IPython", "pandas", "scipy",
    ],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="TADF-Workbench",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,           # GUI app: no console window. Set True to see tracebacks while debugging.
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=None,               # set to "path/to/app.ico" if you add an application icon
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="TADF-Workbench",
)
