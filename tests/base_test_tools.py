# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
#
# This software is a computer program whose purpose is to analyze force curves
# recorded with an atomic force microscope.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use, modify
# and/ or redistribute the software under the terms of the CeCILL license as
# circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading, using, modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean that it is complicated to manipulate, and that also
# therefore means that it is reserved for developers and experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and, more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

"""Some tools for the tests."""

import os
from shutil import copyfileobj
from urllib.request import urlopen

from base_test_files import FOLDER


def dlfile(filename):
    """Downloads the named test file, using urlopen."""

    # Set the source
    url = ("https://bitbucket.org/cmip/pyaf-files/raw/master/tests/" +
           filename)
    # url = ("https://bitbucket.org/cmip/pyaf-files/src/master/tests/" +
    #        filename)

    # Set the destination
    dest = os.path.join(FOLDER, filename)

    # Do not redowmload existing file
    if os.path.isfile(dest):
        print("File " + filename + " already downloaded")
        return

    # Open the url
    with urlopen(url) as in_, open(dest, 'wb') as out:
        print("Downloading " + url)
        copyfileobj(in_, out)

def with_pyplot_interaction(func):
    """Decorator: run the function with pyplot interaction.

    Activate pyplot interaction.
    Needed so that pyplot.show() does not block the system in tests.
    This is not required in pyAF and should not be used because it
    slows down the system.
    """
    def f(*args, **kwargs):
        import matplotlib
        interactive = matplotlib.is_interactive()
        if not interactive:
            matplotlib.interactive(True)
        func(*args, **kwargs)
        if not interactive:
            matplotlib.interactive(False)
    return f
