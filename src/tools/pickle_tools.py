# Copyright Sébastien Janel (2024- ) sebastien.janel@cnrs.fr
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


"""Reading of the pickled data of the saved projects."""

import io
import pickle

# Package of the code: pyAF.src, or src in some scripts
PACKAGE = __name__.rsplit(".tools.", 1)[0]


class ProjectUnpickler(pickle.Unpickler):
    """The classes are saved with their module name, which depends on how
    pyAF was started and on its version: src.experiment or
    pyAF.src.experiment, and experiment for the oldest files. They are looked
    for in the package of the running code: the applications built with
    PyInstaller only have pyAF.src."""

    def find_class(self, module, name):
        for old in ("pyAF.src", "src"):
            if module == old or module.startswith(old + "."):
                module = PACKAGE + module[len(old):]
                break
        else:
            if module == "experiment":
                module = PACKAGE + ".experiment"
        return super().find_class(module, name)


def loads(data, encoding="ASCII"):
    """As pickle.loads, for the data of the saved projects."""
    return ProjectUnpickler(io.BytesIO(data), encoding=encoding).load()
