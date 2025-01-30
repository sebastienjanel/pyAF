#Development notes


##pyAF4
*****************************
This is version of pyAF4 derives from the version ported to python 3 and VTK 5 by PhD.Antoine Dujardin.
It is updated to maintain its compatibilities with new libraries.

Initial Code an main Work by Popoff during his PhD.
Dujardin's Code at https://bitbucket.org/toinedujardin/pyaf/src/refactor/
* Hertz-Fit implemented using [numpy.optimize.curve_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html).
* Added multiprocessing support for Mac OS Catalina. Forkserver method is required for MacOS Catalina, spawn did not seem to work. 
See [bpo-33725](https://bugs.python.org/issue33725) for more details.