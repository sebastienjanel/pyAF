# PyInstaller spec file for building a standalone macOS application
# -*- mode: python ; coding: utf-8 -*-

import os
import sys
import glob

block_cipher = None

def list_files(directory, target_folder):
    """ List all files in a directory and map them to the target folder in dist. """
    return [(os.path.join(directory, file), os.path.join(target_folder, file))
            for file in os.listdir(directory) if os.path.isfile(os.path.join(directory, file))]

# Collect theme styles
theme_files = list_files("./src/gui_styles/stylesheets/dark", "gui_styles/stylesheets/")
theme_files += list_files("./src/gui_styles/stylesheets/light", "gui_styles/stylesheets/")

# Define Analysis
a = Analysis(
    ['main.py'],
    pathex=['.'],
    binaries=[],
    datas=[
        # Images
        (os.path.join(os.getcwd(), "src/images/fix.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/01.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/02.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/03.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/04.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/05.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/06.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/07.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/triangle_down.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/triangle_right.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/icon.icns"), "images/"),
        (os.path.join(os.getcwd(), "src/images/icon.png"), "images/"),
        (os.path.join(os.getcwd(), "src/images/icon.ico"), "images/"),
        # Fonts
        (os.path.join(os.getcwd(), "src/fonts/Lucida_Grande.ttf"), "fonts/"),

        # Stylesheets
        (os.path.join(os.getcwd(), "src/gui_styles/plot_styles/dark_plot.mplstyle"), "gui_styles/plot_styles/"),
        (os.path.join(os.getcwd(), "src/gui_styles/stylesheets/dark.qss"), "gui_styles/stylesheets/"),
        (os.path.join(os.getcwd(), "src/gui_styles/stylesheets/light.qss"), "gui_styles/stylesheets/"),

        # ITK/VTK
        (os.path.join(os.getcwd(), ".venv/lib/python3.13/site-packages/itk"), "itk/"),
        (os.path.join(os.getcwd(), ".venv/lib/python3.13/site-packages/vtkmodules"), "itk/"),
    ] + theme_files,
    hiddenimports=[],
    hookspath=[],
    runtime_hooks=[],
    excludes=[],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False
)

# Build PyInstaller components
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='main',
    debug=True,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='main'
)

app = BUNDLE(
    coll,
    name='pyAF4.app',
    icon='./src/images/icon.icns',
    bundle_identifier=None,
    info_plist={'NSHighResolutionCapable': 'True'}
)

# Create DMG (if needed)
# subprocess.run("create-dmg --overwrite './dist/pyAF4.app' ./dist", shell=True, check=True)
