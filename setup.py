# setup.py is used by py2app to build a standalone application
# but version 0.28 fails in bundling ITK correctly
# pyintaller is the alternative

from setuptools import setup

APP = ['src/pyAF.py']  # Entry point of your application
DATA_FILES = []

OPTIONS = {
    'includes': ['itk'],  # Add ITK to ensure it’s bundled
    'packages': [
        'PIL', 'PyQt5', 'h5py', 'itk', 'lmfit',
        'matplotlib', 'numpy', 'packaging', 'pandas', 'psutil',
        'qtpy', 'scipy', 'seaborn', 'tables',
        'vtk'
    ],
    'excludes': ['setuptools', 'py2app', 'pytest', 'distutils', 'pip'],
    'resources': ['src'],  # Add additional resources if needed
    'plist': {
        'CFBundleName': 'pyAF',
        'CFBundleDisplayName': 'pyAF',
        'LSArchitecturePriority': ['arm64'],
        'CFBundleIdentifier': 'com.CMPI.pyAF',
        'CFBundleVersion': '4.0.0',
        'CFBundleShortVersionString': '4.0.0',
        'NSHumanReadableCopyright': '© 2025 CMPI Lab. All rights reserved.',
    },
}

setup(
    app=APP,
    data_files=DATA_FILES,
    options={'py2app': OPTIONS},
    install_requires=['py2app'],
)