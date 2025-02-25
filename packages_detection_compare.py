import sys
import stdlib_list

# Get the list of built-in modules for your Python version
stdlib_modules = set(stdlib_list.stdlib_list("3.12"))  # Use the closest supported version

# List of imported packages (replace with your actual output)
imported_packages = {
    'PIL', 'PyQt5', 'altgraph', 'anyio', 'asteval', 'blosc2', 'certifi', 'contourpy',
    'cpuinfo', 'cycler', 'dill', 'fontTools', 'h5py', 'h11', 'httpx', 'idna', 'itk',
    'itkConfig', 'kiwisolver', 'lmfit', 'macholib', 'markdown', 'matplotlib', 'modulegraph',
    'mpmath', 'msgpack', 'multiprocessing', 'ndindex', 'numexpr', 'numpy', 'packaging',
    'pandas', 'psutil', 'py2app', 'pyparsing', 'pytz', 'qtpy', 'qtwidgets', 'scipy',
    'seaborn', 'setuptools', 'six', 'sniffio', 'tables', 'typing_extensions',
    'uncertainties', 'vtk', 'vtkmodules'
}

# Remove standard library modules
third_party_packages = imported_packages - stdlib_modules

print("Third-party packages needed:", sorted(third_party_packages))
