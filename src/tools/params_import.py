# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
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

"""Load parameters: the settings of a file of an exported parameters file
(see tools.params_export) given to the files, to analyse new files with the
same parameters. The results are not changed: the files are computed again
with the loaded settings."""

import json

from .params_export import SETTINGS

# Parts of the settings, in the order of the dialog: (key, name)
PARTS = (("elasticity", "Elasticity"),
         ("viscoelasticity", "Viscoelasticity"),
         ("work_and_rupture_force", "Detachment work and force"),
         ("events", "Events"),
         ("loading_rates", "Loading rates"),
         ("breakthrough", "Breakthroughs"),
         ("curves", "Curves (smoothing, tilt)"),
         ("calibration", "Calibration (spring constant, deflection "
                         "sensitivity)"))


# Parts of the computations (the curves and calibration are always given)
COMPUTATIONS = ("elasticity", "viscoelasticity", "work_and_rupture_force",
                "events", "loading_rates", "breakthrough")


class ParametersError(Exception):
    """The file is not a parameters file of pyAF."""


def read(path):
    """Settings of the files of a parameters file: list of (name,
    settings by part)."""
    try:
        with open(path, encoding="utf-8") as f:
            params = json.load(f)
    except (OSError, ValueError) as error:
        raise ParametersError("The file could not be read (" + str(error) +
                              ").") from error
    if not isinstance(params, dict) or "pyaf" not in params or \
            not isinstance(params.get("files"), list):
        raise ParametersError("This is not a parameters file of pyAF.")
    files = [(f.get("name", "?"), f["settings"]) for f in params["files"]
             if isinstance(f, dict) and isinstance(f.get("settings"), dict)]
    if not files:
        raise ParametersError(
            "This parameters file has no settings to load (exported by an "
            "older version of pyAF).")
    # The files with a computation (the others have no settings of
    # analysis, only their calibration and curves)
    computed = [(name, settings) for name, settings in files
                if any(settings.get(part) for part in COMPUTATIONS)]
    if not computed:
        raise ParametersError("No file of this export was computed: no "
                              "settings of analysis to load.")
    return computed


def available_parts(settings):
    """Parts with settings (computed in the exported file)."""
    return [key for key, _ in PARTS if settings.get(key)]


def distinct(files):
    """The files with different settings (one per set of settings)."""
    result = []
    seen = []
    for name, settings in files:
        key = json.dumps(settings, sort_keys=True)
        if key not in seen:
            seen.append(key)
            result.append((name, settings))
    return result


def apply(data, settings, parts):
    """Gives the settings of the parts to a file. Returns the changed
    settings: list of (name, old value, new value)."""
    changes = []
    for part in parts:
        values = settings.get(part) or {}
        for name in SETTINGS[part]:
            if name not in values:
                continue
            new = values[name]
            if part == "calibration" and name == "spring_constant":
                # Stored in N/m, the property gives N/nm
                old = data.spring_constant * 1e9
            else:
                old = getattr(data, name)
            if isinstance(old, bool) or isinstance(new, bool):
                new = bool(new)
            elif isinstance(old, int) and isinstance(new, float) and \
                    new.is_integer():
                new = int(new)
            if old != new:
                changes.append((name, old, new))
                setattr(data, name, new)
    return changes
