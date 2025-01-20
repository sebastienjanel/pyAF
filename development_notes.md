#Development notes


##pyAF3 v.2.0.0
*****************************
This is version of pyAF3 deribes from the version ported to python 3 and VTK 5 by PhD.Antoine Dujardin.

Code: https://bitbucket.org/toinedujardin/pyaf/src/refactor/

This new major version update includes the following features:

* Hertz-Fit implemented using [numpy.optimize.curve_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html).
* Code packaging implemented using [PyInstaller](https://pyinstaller.readthedocs.io/en/stable/). Due to issues for Mac OS X
 support in the latest stable version of PyInstaller (3.6), the development version is used.
* Added multiprocessing support for Mac OS Catalina. Forkserver method is required for MacOS Catalina, spawn did not seem to work. 
See [bpo-33725](https://bugs.python.org/issue33725) for more details.