#Development notes

##pyAF
*****************************
Version 1: initial Code and main work made by Popoff during his PhD (https://theses.fr/api/v1/document/2014LIL10220).

Version 2: ported to python 3 and VTK 5 by PhD. Antoine Dujardin.

Version 3: Hertz-Fit implemented by Javier Lopez-Alonso using [numpy.optimize.curve_fit](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html).
Added multiprocessing support for Mac OS Catalina. Forkserver method is required for MacOS Catalina, spawn did not seem to work. 
See [bpo-33725](https://bugs.python.org/issue33725) for more details.

Version 4: updated by Sébastien Janel to maintain its compatibilities with new libraries.