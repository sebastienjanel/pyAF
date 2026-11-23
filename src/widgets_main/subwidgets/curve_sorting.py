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

"""Curve sorting of the data tab: the curve displayed marked as good or
bad (in any file, map or single curve), and the other curves of all the
files sorted from the marks (see tools/curve_sorting)."""

import logging

import numpy
from qtpy import QtCore, QtWidgets

from ... import shared
from ... import widgets_list
from ...tools import curve_sorting
from ...tools.gui_tools import PYAFButton
from ...widgets.progressbar import Progressbar


class CurveSortingWidget(QtWidgets.QGroupBox):
    """Good and Bad (marks of the curve displayed), Sort all (the curves
    of all the files sorted, the bad ones discarded) and Undo."""

    def __init__(self, parent):
        super().__init__("Curve sorting")
        self.parent = parent
        self.logger = logging.getLogger()
        # Marks, for the session: (unique id of the file, i, j) -> True if
        # bad, False if good
        self.marks = {}
        # Curves discarded by Sort all: unique id of the file -> [(i, j)]
        self.discarded = {}
        # Measures of the curves of the files: unique id -> measures
        self.measures = {}

        self.BT_good = PYAFButton(self, "good", "Good", size=72)
        self.BT_good.setToolTip(
            "Mark the curve displayed as good (G), again to remove the mark")
        self.BT_bad = PYAFButton(self, "bad", "Bad", size=72)
        self.BT_bad.setToolTip(
            "Mark the curve displayed as bad (B), again to remove the mark")
        for button in (self.BT_good, self.BT_bad):
            button.setCheckable(True)
        self.BT_complete = PYAFButton(self, "complete", "Sort all", size=0)
        self.BT_undo = PYAFButton(self, "undo", "Undo", size=0)
        self.BT_undo.setToolTip(
            "Keep again the curves discarded by Sort all")

        # Compact: small gaps between the buttons
        layout = QtWidgets.QHBoxLayout()
        layout.setContentsMargins(6, 2, 6, 4)
        layout.setSpacing(4)
        layout.addStretch(1)
        layout.addWidget(self.BT_good)
        layout.addWidget(self.BT_bad)
        layout.addSpacing(8)
        layout.addWidget(self.BT_complete)
        layout.addWidget(self.BT_undo)
        layout.addStretch(1)
        self.setLayout(layout)
        # As wide as the calibration above
        self.setFixedWidth(350)

        self.update_widget()

    def current(self):
        """File and position of the curve displayed."""
        data = shared.exp.current_data
        if data.is_single:
            return data, 0, 0
        return data, int(shared.exp.meshgrid_click_xpos), \
            int(shared.exp.meshgrid_click_ypos)

    def marks_of(self, data):
        """Marks of a file: [(i, j, True if bad)]."""
        return [(i, j, bad) for (uid, i, j), bad in self.marks.items()
                if uid == data.unique_id]

    def update_widget(self):
        """Counts of the marks, and the mark of the curve displayed."""
        self.forget_removed_files()
        data, i, j = self.current()
        mark = self.marks.get((data.unique_id, i, j))
        bad = sum(self.marks.values())
        good = len(self.marks) - bad
        self.BT_good.setText("Good (%d)" % good if good else "Good")
        self.BT_bad.setText("Bad (%d)" % bad if bad else "Bad")
        self.BT_good.setChecked(mark is False)
        self.BT_bad.setChecked(mark is True)
        tip = "Sort the curves of all the files from the curves marked\n" \
            "Good and Bad (in any file): the bad ones are discarded.\n" \
            "Mark 5 to 10 curves of each kind (the good ones, and each\n" \
            "kind of bad curve), then Sort all; mark the curves sorted\n" \
            "wrongly, and Sort all again."
        if not curve_sorting.available():
            tip = "Needs scikit-learn (not installed)"
        self.BT_complete.setToolTip(tip)
        self.BT_complete.setEnabled(curve_sorting.available() and
                                    bool(self.marks))
        self.BT_undo.setEnabled(any(self.discarded.values()))

    def forget_removed_files(self):
        """No marks, discarded curves or measures of the removed files."""
        ids = {data.unique_id for data in shared.exp.list}
        self.marks = {key: bad for key, bad in self.marks.items()
                      if key[0] in ids}
        for store in (self.discarded, self.measures):
            for uid in [uid for uid in store if uid not in ids]:
                del store[uid]

    def button_clicked(self, name):
        """Called when a button is clicked."""
        if name in ("good", "bad"):
            self.mark(name == "bad")
        elif name == "complete":
            self.complete()
        elif name == "undo":
            self.undo()
            self.refresh()
            self.update_widget()

    def mark(self, bad):
        """The curve displayed marked as bad (or good), or its mark removed
        if it has this one."""
        data, i, j = self.current()
        if data.is_outside(i, j):
            return
        key = (data.unique_id, i, j)
        if self.marks.get(key) is bad:
            del self.marks[key]
        else:
            self.marks[key] = bad
        self.logger.debug("Curve %s marked %s", key, self.marks.get(key))
        self.update_widget()
        if not data.is_single:
            canvas = widgets_list.widget_data.meshgrid.canvas
            if hasattr(canvas, "update_blit"):
                canvas.update_blit("all")

    def complete(self):
        """The curves of all the files sorted from the marks, the bad ones
        discarded (and the ones discarded by the last Sort all kept
        again)."""
        self.forget_removed_files()
        files = list(shared.exp.list)
        index = {data.unique_id: number for number, data in enumerate(files)}
        marks = [(index[uid], i, j, bad)
                 for (uid, i, j), bad in self.marks.items()]
        missing = [data for data in files if data.unique_id not in
                   self.measures]
        if missing:
            Progressbar("Curve sorting")
            bar = widgets_list.widget_progressbar
            bar.set_label("Reading the curves")
            bar.set_range(0, len(missing))
            try:
                for data in missing:
                    self.measures[data.unique_id] = \
                        curve_sorting.file_measures(data)
                    bar.update()
            finally:
                if widgets_list.widget_progressbar is not None:
                    widgets_list.widget_progressbar.close()
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            all_features = curve_sorting.features(
                [self.measures[data.unique_id] for data in files],
                [data.get_fit_quality(0) if data.has_fit_quality() else None
                 for data in files])
            probabilities = curve_sorting.bad_probability(all_features,
                                                          marks)
        except curve_sorting.CurveSortingError as error:
            QtWidgets.QApplication.restoreOverrideCursor()
            QtWidgets.QMessageBox.information(self, "Curve sorting",
                                              str(error))
            return
        self.undo()
        total = 0
        for data, proba in zip(files, probabilities):
            bad = numpy.nan_to_num(proba, nan=0.0) >= 0.5
            for (uid, i, j), marked_bad in self.marks.items():
                if uid == data.unique_id:
                    bad[i, j] = marked_bad
            curves = []
            for i, j in numpy.argwhere(bad):
                if data.discarded_curves[i][j] == 0 and \
                        not data.is_outside(i, j):
                    data.discarded_curves[i][j] = 1
                    curves.append((int(i), int(j)))
            self.discarded[data.unique_id] = curves
            total += len(curves)
        QtWidgets.QApplication.restoreOverrideCursor()
        self.logger.debug("Curve sorting: %s curves discarded", total)
        self.refresh()
        self.update_widget()
        number = sum(1 for curves in self.discarded.values() if curves)
        QtWidgets.QMessageBox.information(
            self, "Curve sorting",
            "%d curve%s discarded in %d file%s (Undo keeps them again)."
            % (total, "s" if total != 1 else "", number,
               "s" if number != 1 else ""))

    def undo(self):
        """The curves discarded by the last Sort all kept again (not the
        ones kept or discarded otherwise since)."""
        files = {data.unique_id: data for data in shared.exp.list}
        for uid, curves in self.discarded.items():
            data = files.get(uid)
            if data is None:
                continue
            for i, j in curves:
                data.discarded_curves[i][j] = 0
        self.discarded = {}

    @staticmethod
    def refresh():
        """Maps, curve and results with the curves discarded."""
        data = shared.exp.current_data
        if not data.is_single:
            canvas = widgets_list.widget_data.meshgrid.canvas
            if hasattr(canvas, "update_blit"):
                canvas.update_blit("all")
        widgets_list.widget_data.curve.update_plot()
        results = widgets_list.widget_results
        if results is not None and not data.is_single:
            canvas = getattr(results.MPL_canvas1, "canvas", None)
            if hasattr(canvas, "update_blit"):
                canvas.update_blit("discarded")
        if widgets_list.widget_curve_mod is not None:
            widgets_list.widget_curve_mod.update_widget()
        # The values of the results without the curves discarded
        main = widgets_list.widget_main
        if main is not None:
            if main.tabs.currentIndex() in (3, 4):
                widgets_list.widget_results_single.schedule_refresh()
            else:
                main.update_results = True
