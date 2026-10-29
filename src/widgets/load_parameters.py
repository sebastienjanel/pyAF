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

"""Dialog of Load parameters (Compute tab): which settings of an exported
parameters file are given to the current file (or to all the files with
"Apply to all"), see tools.params_import."""

import os

from qtpy import QtCore, QtWidgets

from .. import shared
from ..tools import gui_tools, misc_tools, params_import
from ..tools.gui_tools import PYAFComboBox, PYAFToggle


def load_parameters(parent):
    """Asks for a parameters file and the settings to load, and gives them
    to the files. Returns True if settings were loaded."""
    settings = QtCore.QSettings()
    path, _ = QtWidgets.QFileDialog.getOpenFileName(
        parent, "Load parameters", misc_tools.get_user_path(settings),
        "pyAF parameters (*.json)")
    if not path:
        return False
    misc_tools.set_user_path(settings, os.path.dirname(path))
    try:
        files = params_import.read(path)
    except params_import.ParametersError as error:
        QtWidgets.QMessageBox.warning(parent, "Load parameters", str(error))
        return False

    dialog = LoadParametersDialog(parent, path, files)
    if dialog.exec_() != QtWidgets.QDialog.Accepted:
        return False
    file_settings, parts = dialog.choice()
    if not parts:
        return False

    if shared.exp.apply_to_all_compute:
        targets = list(shared.exp.list)
    else:
        targets = [shared.exp.current_data]
    changes = {}
    for data in targets:
        for name, old, new in params_import.apply(data, file_settings,
                                                  parts):
            changes.setdefault(name, (old, new))

    target = "all the files" if len(targets) > 1 else \
        os.path.basename(str(targets[0].filename))
    if changes:
        lines = ["%s: %s → %s" % (name, old, new)
                 for name, (old, new) in sorted(changes.items())]
        gui_tools.show_success(
            parent, "Parameters loaded for " + target + ": " +
            str(len(changes)) + " setting" + ("s" if len(changes) > 1
                                              else "") +
            " changed.\nCompute again to apply them.\n\n" + "\n".join(lines))
    else:
        gui_tools.show_success(
            parent, "Parameters loaded for " + target + ": the settings "
            "were already the same.")
    return True


class LoadParametersDialog(QtWidgets.QDialog):
    """Choice of the file of the parameters (if they differ) and of the
    parts of the settings to load."""

    def __init__(self, parent, path, files):
        super().__init__(parent)
        self.setWindowTitle("Load parameters")
        self.files = params_import.distinct(files)

        VL = QtWidgets.QVBoxLayout()
        VL.addWidget(QtWidgets.QLabel(
            "Parameters of <b>" + os.path.basename(path) + "</b>"))

        # Files with different settings: choice of one
        self.list_files = None
        if len(self.files) > 1:
            HL = QtWidgets.QHBoxLayout()
            HL.addWidget(QtWidgets.QLabel("Settings of"))
            self.list_files = PYAFComboBox(self, "file")
            for name, _ in self.files:
                self.list_files.addItem(name)
            self.list_files.fit_to_contents()
            HL.addWidget(self.list_files)
            HL.addStretch(1)
            VL.addLayout(HL)
            note = QtWidgets.QLabel(
                "The files of this export were computed with different "
                "settings.")
            note.setWordWrap(True)
            VL.addWidget(note)

        # Parts of the settings (the calibration belongs to the cantilever:
        # not loaded by default)
        self.toggles = {}
        GL = QtWidgets.QGridLayout()
        width = max(QtWidgets.QLabel(name).sizeHint().width()
                    for _, name in params_import.PARTS)
        for row, (key, name) in enumerate(params_import.PARTS):
            toggle = PYAFToggle(self, key, name)
            toggle.label.setFixedWidth(width + 8)
            self.toggles[key] = toggle
            GL.addWidget(toggle, row, 0, QtCore.Qt.AlignLeft)
        box = QtWidgets.QGroupBox("Settings to load")
        box.setLayout(GL)
        VL.addWidget(box)

        target = "all the files (Apply to all)" \
            if shared.exp.apply_to_all_compute else "the current file (" + \
            os.path.basename(str(shared.exp.current_data.filename)) + ")"
        info = QtWidgets.QLabel(
            "Given to " + target + ". The results do not change until they "
            "are computed again.")
        info.setWordWrap(True)
        VL.addWidget(info)

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        VL.addWidget(buttons)
        self.setLayout(VL)
        self.update_toggles()

    def current_settings(self):
        """Settings of the chosen file."""
        index = self.list_files.currentIndex() if self.list_files else 0
        return self.files[index][1]

    def update_toggles(self):
        """Parts of the chosen file: the computed ones, on (but the
        calibration), the others off and disabled."""
        available = params_import.available_parts(self.current_settings())
        for key, toggle in self.toggles.items():
            toggle.setEnabled(key in available)
            toggle.setChecked(key in available and key != "calibration")

    def list_updated(self, _):
        """Another file was chosen."""
        self.update_toggles()

    def toggle_clicked(self, _):
        """A part was (un)chosen."""

    def choice(self):
        """Settings of the chosen file and chosen parts."""
        parts = [key for key, toggle in self.toggles.items()
                 if toggle.isEnabled() and toggle.isChecked()]
        return self.current_settings(), parts
