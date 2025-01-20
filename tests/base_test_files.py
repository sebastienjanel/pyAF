# Copyright Antoine Dujardin (2016-2018) toine.dujardin@gmail.com
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

import os.path

from pyAF.src.tools import misc_tools


def _def_file(filename, file_type, nbrcurves):
    """Return the dict description of the file.

    The following arguments must be provided:
        - filename: name of the file;
        - file_type: see below;
        - nbrcurves: number of curves.
    """
    return {
        "version": None,
        "checked": True,
        "error": "",
        "file_type": file_type,
        "path": os.path.join(FOLDER, filename),
        "nbrcurves": 1,
        "enabled": True,
        "filename": filename
    }

FOLDER = os.path.join(
    misc_tools.get_app_path(), os.path.pardir, "tests", "files")

# file_type
_NS_SF_T = "Nanoscope (Single File)"
_NS_FV_T = "Nanoscope (Force Volume)"
_JPK_SF_T = "JPK (Single File)"
_JPK_FM_T = "JPK (Force Map)"
_JPK_QI_T = "JPK (QI)"

HELA = _def_file("hela.003", _NS_FV_T, 256)
NS_SINGLE = _def_file("single_curve.001", _NS_SF_T, 1)
JPK_SINGLE = _def_file("43_curve1.jpk-force", _JPK_SF_T, 1)
PTK2 = _def_file("ptk2.000", _NS_FV_T, 256)
MAP16 = _def_file("map16ums.jpk-force-map", _JPK_FM_T, 1024)
NONSQR1 = _def_file("nonsquare1.jpk-qi-data", _JPK_QI_T, 1827)
NONSQR2 = _def_file("nonsquare2.jpk-qi-data", _JPK_QI_T, 1827)
NS130226001 = _def_file("130226.001", _NS_FV_T, 256)
NS130226002 = _def_file("130226.002", _NS_FV_T, 256)
NS130226003 = _def_file("130226.003", _NS_FV_T, 256)
JPK_S_SPE = _def_file("force-save-act4-2014.06.16-17.13.21.jpk-force",
                      _JPK_SF_T, 1)
JPK_S_RAMP_SPE = _def_file("force-save-2014.06.20-15.21.59.jpk-force",
                           _JPK_SF_T, 1)
JPK_OLD_MAP = _def_file("map-data-2012SJ.jpk-force-map", _JPK_FM_T, 1)
JPK_CORRUPTED = _def_file("corrupted_pixels.jpk-force-map", _JPK_FM_T, 256)
