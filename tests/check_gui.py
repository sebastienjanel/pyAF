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


"""Checks the Data, Compute and Results tabs of pyAF with your own files.

Loads the files, computes the elasticity, the work and rupture forces, the
unbinding events and the loading rates on all of them, then uses every
control of the three tabs (buttons, checkboxes, radio buttons, lists and
input fields), for each file, before and after the computations.

Reports the errors (with the action which triggered them), the Qt messages
and the Python warnings coming from pyAF. The exit code is 1 if there are
errors.

Usage (from the pyAF folder):
    .venv/bin/python tests/check_gui.py file1 [file2 ...]

Use files of different types (e.g. a JPK force map, a JPK QI map and a
Bruker file) to check them together. By default no window is displayed
(QT_QPA_PLATFORM=offscreen); the 3D view is not checked, as it needs a
display. Your pyAF preferences are not modified.
"""

import collections
import faulthandler
import os
import sys
import tempfile
import traceback
import warnings

os.environ.setdefault("QT_API", "pyside6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(ROOT))
sys.path.insert(0, ROOT)

# Buttons not clicked: file management, and the 3D view (needs a display)
SKIPPED_BUTTONS = {"+", "-", "Rename", "Copy", "Display in 3D"}


class Report:
    """Collects the errors, dialogs, Qt messages and warnings."""

    def __init__(self):
        self.action = "start"
        self.errors = []
        self.dialogs = collections.Counter()
        self.qt_messages = collections.Counter()
        self.warnings = collections.Counter()

    def exception(self, exc_type, value, tb):
        text = "".join(traceback.format_exception(exc_type, value, tb))
        self.errors.append((self.action, text))

    def warning(self, message, category, filename, lineno, file=None,
                line=None):
        if os.sep + "src" + os.sep in filename and ROOT in filename:
            where = os.path.relpath(filename, ROOT) + ":" + str(lineno)
            self.warnings[category.__name__ + ": " + str(message)[:100] +
                          " (" + where + ")"] += 1

    def qt_message(self, mode, context, message):
        if "propagateSizeHints" in message:
            # Message of the offscreen platform, not of pyAF
            return
        self.qt_messages[message[:140]] += 1

    def dialog(self, text):
        self.dialogs["[" + self.action[:60] + "] " + text[:90]] += 1

    def print(self):
        print("\n=== Errors:", len(self.errors))
        by_error = collections.OrderedDict()
        for action, text in self.errors:
            last_line = text.strip().splitlines()[-1]
            by_error.setdefault(last_line, []).append((action, text))
        for last_line, occurrences in by_error.items():
            print("\n--- " + last_line + " (x" + str(len(occurrences)) + ")")
            actions = sorted(set(action for action, _ in occurrences))
            for action in actions[:10]:
                print("    after: " + action)
            print(occurrences[0][1])

        print("\n=== Dialogs:", sum(self.dialogs.values()))
        for text, count in self.dialogs.most_common(30):
            print("  x" + str(count), text)

        print("\n=== Qt messages:", len(self.qt_messages))
        for text, count in self.qt_messages.most_common(20):
            print("  x" + str(count), text)

        print("\n=== Python warnings:", len(self.warnings))
        for text, count in self.warnings.most_common(20):
            print("  x" + str(count), text)


def main(paths):
    from qtpy import QtCore, QtWidgets
    from qtpy.QtTest import QTest

    report = Report()
    faulthandler.enable()
    warnings.simplefilter("always")
    warnings.showwarning = report.warning
    QtCore.qInstallMessageHandler(report.qt_message)

    app = QtWidgets.QApplication(["pyAF"])

    # Settings separated from the user's ones
    QtCore.QCoreApplication.setOrganizationName("pyAF-check-gui")
    QtCore.QCoreApplication.setApplicationName("pyAF-check-gui")
    temp_folder = tempfile.mkdtemp(prefix="pyaf_check_gui_")
    settings = QtCore.QSettings()
    settings.clear()
    settings.setValue("tempPath", temp_folder + os.sep)

    # Never block on a dialog: answer No / Cancel
    def message_box_exec(box):
        report.dialog("QMessageBox: " + box.text())
        return QtWidgets.QMessageBox.No

    QtWidgets.QMessageBox.exec = message_box_exec
    QtWidgets.QMessageBox.exec_ = message_box_exec
    for name in ("warning", "information", "critical", "question"):
        def static_box(*args, name=name, **kwargs):
            report.dialog(name + ": " + (str(args[2]) if len(args) > 2
                                         else ""))
            return QtWidgets.QMessageBox.No
        setattr(QtWidgets.QMessageBox, name, staticmethod(static_box))
    QtWidgets.QInputDialog.getText = staticmethod(lambda *a, **k: ("", False))
    QtWidgets.QFileDialog.getOpenFileNames = \
        staticmethod(lambda *a, **k: ([], ""))
    QtWidgets.QFileDialog.getSaveFileName = \
        staticmethod(lambda *a, **k: ("", ""))
    QtWidgets.QFileDialog.getExistingDirectory = \
        staticmethod(lambda *a, **k: "")

    def close_modal_dialogs():
        for widget in QtWidgets.QApplication.topLevelWidgets():
            if isinstance(widget, QtWidgets.QDialog) and \
                    widget.isVisible() and widget.isModal():
                report.dialog(type(widget).__name__ + ": " +
                              widget.windowTitle())
                widget.reject()

    timer = QtCore.QTimer()
    timer.timeout.connect(close_modal_dialogs)
    timer.start(150)

    from src import consts, shared, widgets_list, pyaf
    from src.load_and_save.pre_check import CheckFiles
    from src.gui_styles.theme_handler import theme_handler

    consts.UNIT_TESTING = True
    sys.excepthook = report.exception
    theme_handler("Dark")

    def act(name, function, wait=60):
        report.action = name
        try:
            function()
        except Exception:
            report.exception(*sys.exc_info())
        QTest.qWait(wait)

    main_window = pyaf.MainWindow(CheckFiles(paths).infos)
    main_window.show()
    QTest.qWait(1200)

    import qtpy
    print("Qt:", qtpy.API_NAME, qtpy.QT_VERSION)
    for data in shared.exp.list:
        print("File:", data.filename, "(" + data.file_type + ")")

    def exercise(widget, label):
        """Use every enabled control below widget. Returns the number of
        actions."""
        actions = 0
        for button in widget.findChildren(QtWidgets.QAbstractButton):
            if not (button.isVisible() and button.isEnabled()):
                continue
            text = button.text().strip() or button.objectName()
            if text in SKIPPED_BUTTONS:
                continue
            if button.isCheckable():
                act(label + ": toggle '" + text + "'", button.click)
                actions += 1
                if isinstance(button, QtWidgets.QCheckBox) and \
                        button.isVisible() and button.isEnabled():
                    act(label + ": toggle back '" + text + "'", button.click)
                    actions += 1
            else:
                act(label + ": click '" + text + "'", button.click, 120)
                actions += 1
                # Close the windows opened by the button
                for window in QtWidgets.QApplication.topLevelWidgets():
                    if window is not main_window and window.isVisible() \
                            and window.isWindow() and \
                            not isinstance(window, QtWidgets.QMenu):
                        act(label + ": close the window opened by '" +
                            text + "'", window.close)

        for combo in widget.findChildren(QtWidgets.QComboBox):
            if not (combo.isVisible() and combo.isEnabled()):
                continue
            current = combo.currentIndex()
            for index in list(range(combo.count())) + [current]:
                if combo.isVisible() and combo.isEnabled():
                    act(label + ": list '" + combo.objectName() + "' -> '" +
                        combo.itemText(index) + "'",
                        lambda c=combo, i=index: c.setCurrentIndex(i))
                    actions += 1

        for line_edit in widget.findChildren(QtWidgets.QLineEdit):
            if not (line_edit.isVisible() and line_edit.isEnabled()) or \
                    line_edit.isReadOnly() or \
                    isinstance(line_edit.parent(), QtWidgets.QComboBox):
                continue
            old = line_edit.text()
            # 0 (refused by some inputs), then the original value
            for value in ("0", old):
                def set_value(e=line_edit, v=value):
                    e.setText(v)
                    if e.hasAcceptableInput():
                        e.editingFinished.emit()
                act(label + ": input '" + old + "' -> '" + value + "'",
                    set_value)
                actions += 1
        return actions

    tabs = {0: ("Data", widgets_list.widget_data),
            1: ("Compute", widgets_list.widget_compute),
            2: ("Results", widgets_list.widget_results)}
    compute = widgets_list.widget_compute

    def check_tabs(step):
        for tab_index, (tab_name, tab_widget) in tabs.items():
            for file_index, data in enumerate(shared.exp.list):
                act("select file " + str(file_index),
                    lambda i=file_index:
                    widgets_list.widget_data.list_boxes[i].change_id(), 400)
                act("show the " + tab_name + " tab",
                    lambda i=tab_index: main_window.tabs.setCurrentIndex(i),
                    300)
                label = tab_name + ", " + data.filename + ", " + step
                actions = 0
                if tab_name == "Compute":
                    for sub_tab in range(compute.tabs.count()):
                        act("show compute tab " + str(sub_tab),
                            lambda i=sub_tab: compute.tabs.setCurrentIndex(i),
                            200)
                        actions += exercise(
                            tab_widget, label + ", " +
                            compute.tabs.tabText(sub_tab))
                else:
                    actions += exercise(tab_widget, label)
                print(label + ": " + str(actions) + " actions, " +
                      str(len(report.errors)) + " errors so far",
                      flush=True)

    check_tabs("before computing")

    main_window.tabs.setCurrentIndex(1)
    QTest.qWait(300)
    shared.exp.calc_all = True
    for sub_tab, button in ((0, compute.BT_compute_stiffness),
                            (2, compute.BT_compute_work_and_rupture_force),
                            (3, compute.BT_compute_events)):
        act("compute " + compute.tabs.tabText(sub_tab),
            lambda i=sub_tab, b=button: (compute.tabs.setCurrentIndex(i),
                                         QTest.mouseClick(
                                             b, QtCore.Qt.LeftButton)),
            1500)
    act("compute the loading rates",
        lambda: (compute.tabs.setCurrentIndex(4),
                 compute.button_clicked("get_lr")), 1500)
    for data in shared.exp.list:
        print(data.filename + ": elasticity", data.stiffness_calculated,
              ", work", data.work_and_rupture_force1_calculated,
              ", events", data.events_calculated)

    check_tabs("after computing")

    report.print()
    main_window.close()
    return 1 if report.errors else 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main([os.path.abspath(path) for path in sys.argv[1:]]))
