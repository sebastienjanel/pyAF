from setuptools import setup

APP = ['src/pyAF.py']  # Entry point of your application
DATA_FILES = []

OPTIONS = {
    'includes': ['itk'],  # Add ITK to ensure it’s bundled
    'packages': [
        'PIL', 'PyQt5', 'altgraph', 'anyio', 'asteval', 'blosc2', 'certifi',
        'contourpy', 'cpuinfo', 'cycler', 'dill', 'fontTools', 'h11', 'h5py',
        'httpx', 'idna', 'itk', 'itkConfig', 'kiwisolver', 'lmfit', 'macholib',
        'markdown', 'matplotlib', 'modulegraph', 'mpmath', 'msgpack', 'ndindex',
        'numexpr', 'numpy', 'packaging', 'pandas', 'psutil', 'pyparsing', 'pytz',
        'qtpy', 'qtwidgets', 'scipy', 'seaborn', 'tables', 'typing_extensions',
        'uncertainties', 'vtk', 'vtkmodules'
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