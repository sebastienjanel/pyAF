# PyInstaller spec file: standalone application of pyAF for macOS, Windows
# and Linux (built on each system, PyInstaller does not cross-compile).
#
#     pip install pyinstaller
#     pyinstaller pyAF.spec
#
# Result in dist/: pyAF.app (macOS) or the pyAF folder with the pyAF
# executable (Windows, Linux). See .github/workflows/build-apps.yml.
# -*- mode: python ; coding: utf-8 -*-

import os
import sys

# Version shown by the application (Info.plist on macOS), set by the build
# workflow from the git tag
VERSION = os.environ.get("PYAF_VERSION", "5")

SRC = os.path.join(os.getcwd(), "src")


def resources(folder):
    """Files of a folder of src (not the code), at the same place in the
    bundle (see misc_tools.get_app_path)."""
    files = []
    for root, _, names in os.walk(os.path.join(SRC, folder)):
        if "__pycache__" in root:
            continue
        for name in names:
            if name.endswith((".py", ".pyc")) or name == ".DS_Store":
                continue
            path = os.path.join(root, name)
            files.append((path, os.path.relpath(root, SRC)))
    return files


datas = resources("images") + resources("fonts") + resources("gui_styles")

if sys.platform == "darwin":
    icon = os.path.join(SRC, "images", "icon.icns")
elif sys.platform == "win32":
    icon = os.path.join(SRC, "images", "icon.ico")
else:
    icon = None

a = Analysis(
    ["main.py"],
    # Parent of the repository folder: pyAF imports itself as the pyAF
    # package (see main.py)
    pathex=[os.path.dirname(os.getcwd())],
    binaries=[],
    datas=datas,
    # VTK, PyTables, Qt: hooks of PyInstaller and
    # pyinstaller-hooks-contrib
    hiddenimports=["pyAF.src.pyaf"],
    hookspath=[],
    runtime_hooks=[],
    # Only PySide6 (PyInstaller refuses several Qt bindings), and modules
    # not used by pyAF
    # "vtk" would bring all of VTK (see vtk_api.py)
    excludes=["PyQt5", "PyQt6", "PySide2", "tkinter", "IPython",
              "jupyter", "notebook", "pytest", "vtk", "vtkmodules.all",
              "PySide6.QtNetwork", "PySide6.QtDataVisualization"],
    noarchive=False,
)


def used(entry):
    """False for the parts of Qt not used by pyAF, that nothing needed links
    to (they come with the hooks of PySide6), and for the libraries to take
    from the system."""
    name = entry[0].replace("\\", "/")
    if name in SYSTEM_LIBRARIES:
        return False
    return not any(part in name for part in UNUSED_QT)


# Frameworks (macOS), Qt6 libraries (Windows, Linux) and plugins
UNUSED_QT = ["Qt" + module for module in
             ("Quick", "Qml", "Pdf", "Network", "DataVisualization",
              "VirtualKeyboard")]
UNUSED_QT += ["Qt6" + module for module in
              ("Quick", "Qml", "Pdf", "Network", "DataVisualization",
               "VirtualKeyboard")]
UNUSED_QT += ["plugins/tls/", "plugins/networkinformation/", "tuiotouch",
              "qtvirtualkeyboard", "imageformats/libqpdf",
              "imageformats/qpdf"]
# Linux: the graphics driver comes from the system and needs its C++ library:
# with the older one of the system of the build, the driver can not be loaded
# and the 3D view crashes (seen on Ubuntu 26.04). All the supported systems
# have one at least as recent as the system of the build.
SYSTEM_LIBRARIES = ["libstdc++.so.6", "libgcc_s.so.1"] \
    if sys.platform.startswith("linux") else []
a.binaries = [entry for entry in a.binaries if used(entry)]
a.datas = [entry for entry in a.datas if used(entry)]

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="pyAF",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    # UPX breaks Qt and VTK libraries
    upx=False,
    # Window application (no terminal); on Linux the terminal stays usable
    console=False,
    icon=icon,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=False,
    name="pyAF",
)

if sys.platform == "darwin":
    app = BUNDLE(
        coll,
        name="pyAF.app",
        icon=icon,
        bundle_identifier="fr.cnrs.cmpi.pyaf",
        version=VERSION,
        info_plist={
            "NSHighResolutionCapable": True,
            "CFBundleShortVersionString": VERSION,
            "CFBundleDisplayName": "pyAF",
            "NSRequiresAquaSystemAppearance": False,
        },
    )
