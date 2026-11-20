# Copyright Michka Popoff (2011-2014) michkapopoff@gmail.com
# Copyright Antoine Dujardin (2016-2017) toine.dujardin@gmail.com
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

"""Widget used to manage the regions of interest (ROIs), drawn on its own
map."""

import math

import numpy
from .. import shared
from .. import widgets_list
import logging
from qtpy import QtCore, QtGui, QtWidgets
from ..compute_tools import plant
from ..compute_tools.stiffness import bec
from ..plots.PYAFPlot import PYAFPlot
from ..plots.plot_main import MainPlot
from ..tools import auto_roi
from ..tools import gui_tools
from ..tools import map_values
from ..tools import math_tools
from ..tools import misc_tools
from ..tools import planefit
from ..tools import plant_tools
from ..tools.gui_tools import PYAFButton
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFButtonGroup
from ..tools.PYAFWidget import PYAFWidget
from ..tools.gui_tools import PYAFInput
from ..tools.gui_tools import PYAFToggle
from .. import experiment


class RoiManagerWidget(PYAFWidget):
    """Widget used to manage the regions of interest (ROIs).

    The ROIs are drawn on the map of the window (piezo height as soon as the
    file is loaded, or any computed map), before or after the computations.
    They define the substrate (B.E.C., planefit), the cells (plant cells),
    or sort the data in the results.
    """

    def __init__(self, parent):
        super().__init__(parent, "widget_roi_manager")
        self.setWindowTitle("ROI Manager")

        self.BG_selected = None
        # Drawing: "select" (click: curve), "roi_add" (default) or
        # "roi_erase"
        self.cursor_mode = "roi_add"
        # Tool of the maps while a single curve is shown (Select curve)
        self.map_tool = None
        self.drawing = False
        # File of the map displayed
        self.shown_file = None
        # ROI assistance: pixels marked in each ROI and pixels erased from it
        # (sets of (i, j) by ROI, by file, for the session), features of the
        # file of the map, and the proposal (see auto_update)
        self.auto_marks = {}
        self.auto_erased = {}
        self.auto_features = None
        self.auto_domes = None
        self.auto_result = None

        self.logger = logging.getLogger()
        self.logger.debug("Opening ROI manager")

        self.VL = QtWidgets.QVBoxLayout()
        self.HL = QtWidgets.QHBoxLayout()

        # --- Map -------------------------------------------------------------

        self.list_files = PYAFComboBox(self, "files")
        self.list_files.setToolTip("File of the map (the single curves have "
                                   "no map)")
        self.list_maps = PYAFComboBox(self, "maps")
        self.list_maps.setToolTip(
            "Map to draw on: the piezo height as soon as the file is loaded, "
            "the other maps once computed")
        # File, and the map below it
        GL_lists = QtWidgets.QGridLayout()
        GL_lists.addWidget(QtWidgets.QLabel("File"), 0, 0)
        GL_lists.addWidget(self.list_files, 0, 1)
        GL_lists.addWidget(QtWidgets.QLabel("Map"), 1, 0)
        GL_lists.addWidget(self.list_maps, 1, 1, QtCore.Qt.AlignLeft)
        GL_lists.setColumnStretch(1, 1)

        resolution = widgets_list.widget_main.canvas_resolution
        canvas_size = [520, 460]
        # The map fills the free space of the window (initial size: the
        # size of the window)
        self.W_canvas = QtWidgets.QWidget()
        self.W_canvas.setMinimumSize(canvas_size[0], canvas_size[1])
        self.W_canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                    QtWidgets.QSizePolicy.Expanding)
        sizes = [canvas_size[0] / resolution, canvas_size[1] / resolution,
                 resolution]
        self.MPL_canvas = PYAFPlot(self, "meshgrid_roi", self.W_canvas,
                                   sizes)
        layout_canvas = QtWidgets.QVBoxLayout(self.W_canvas)
        layout_canvas.setContentsMargins(0, 0, 0, 0)
        layout_canvas.addWidget(self.MPL_canvas)
        self.MPL_canvas.mpl_connect("button_press_event", self.canvas_press)
        self.MPL_canvas.mpl_connect("button_release_event",
                                    self.canvas_release)
        self.MPL_canvas.mpl_connect("motion_notify_event", self.canvas_motion)

        # Drawing tools
        self.BG_tool = PYAFButtonGroup(self, "tool")
        self.RB_select = QtWidgets.QRadioButton("Select curve")
        self.RB_select.setToolTip("Click: curve displayed in the main tabs")
        self.RB_draw = QtWidgets.QRadioButton("Draw")
        self.RB_draw.setToolTip("Add the pixels under the pencil to the ROI "
                                "checked in the Draw column")
        self.RB_erase = QtWidgets.QRadioButton("Erase")
        self.RB_erase.setToolTip("Remove the pixels under the pencil from "
                                 "the ROI checked in the Draw column")
        for index, button in enumerate((self.RB_select, self.RB_draw,
                                        self.RB_erase)):
            self.BG_tool.addButton(button, index)
        self.RB_draw.setChecked(True)
        self.TG_auto = PYAFToggle(self, "roi_assistance", "ROI assistance")
        if auto_roi.available():
            self.TG_auto.setToolTip(
                "Sort the other pixels into the ROIs from a few pixels of\n"
                "them (the ROIs drawn, and the pixels marked with Draw\n"
                "while ROI assistance is on): by regions with all the maps\n"
                "computed for the file (topography, elasticity...), or by\n"
                "dome with the domes of the height (see Sort by). Erase\n"
                "takes pixels out of the ROI (by dome, with the region\n"
                "around them up to the grooves). The ROIs change only with\n"
                "Accept.")
        else:
            self.TG_auto.setToolTip("Needs scikit-learn (not installed)")
        self.change_cursor_mode(self.cursor_mode)
        self.label_no_map = QtWidgets.QLabel(
            "Single curve: no map. Its cell is drawn on the map where it was "
            "recorded.")
        self.label_no_map.setWordWrap(True)
        HL_tools = QtWidgets.QHBoxLayout()
        HL_tools.addWidget(self.RB_select)
        HL_tools.addWidget(self.RB_draw)
        HL_tools.addWidget(self.RB_erase)
        HL_tools.addSpacing(12)
        HL_tools.addWidget(self.TG_auto)
        HL_tools.addStretch(1)

        self.VL_map = QtWidgets.QVBoxLayout()
        self.VL_map.addLayout(GL_lists)
        self.VL_map.addWidget(self.W_canvas, 1)
        self.VL_map.addWidget(self.label_no_map)
        self.VL_map.addLayout(HL_tools)
        self.HL.addLayout(self.VL_map)

        self.TW = QtWidgets.QTableWidget(0, 8)
        self.TW.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.TW.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.TW.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        labels = ["#", "Draw", "Show", "Pixels", "Area [µm²]", "Color",
                  "Substrate", "Data"]
        self.TW.setHorizontalHeaderLabels(labels)
        tooltips = ["", "ROI drawn or erased with the Draw and Erase tools",
                    "Display the ROI on the maps", "", "", "",
                    "ROI defined as substrate", "Curves of the ROI "
                    "discarded with Discard data"]
        for col, tooltip in enumerate(tooltips):
            self.TW.horizontalHeaderItem(col).setToolTip(tooltip)
        # At least the width of the titles (larger fonts on Linux and
        # Windows), and the list as wide as its columns
        self.TW.ensurePolished()
        header = self.TW.horizontalHeader()
        header.setMinimumSectionSize(20)
        widths = [30, 45, 45, 60, 85, 100, 75, 80]
        for col, width in enumerate(widths):
            self.TW.setColumnWidth(
                col, max(width, header.sectionSizeHint(col)))
        self.TW.setFixedWidth(header.length() + 2 * self.TW.frameWidth())
        # The last column takes the rest of the width (no empty space on
        # the right)
        header.setStretchLastSection(True)
        self.TW.verticalHeader().hide()
        self.TW.cellClicked.connect(self.tableWidget_clicked)

        self.VL_options = QtWidgets.QVBoxLayout()

        self.input_size = PYAFInput(self, "size", "Pencil size")
        self.BG = PYAFButtonGroup(self, "bg")
        self.radiobutton_square = QtWidgets.QRadioButton("Square")
        self.radiobutton_circle = QtWidgets.QRadioButton("Circle")
        self.radiobutton_dot = QtWidgets.QRadioButton("Dot")
        self.BG.addButton(self.radiobutton_square, 0)
        self.BG.addButton(self.radiobutton_circle, 1)
        self.BG.addButton(self.radiobutton_dot, 2)

        self.button_add = PYAFButton(self, "add", "Add ROI")
        self.button_remove = PYAFButton(self, "remove", "Remove")
        self.button_invert = PYAFButton(self, "invert", "Invert to new")
        self.BT_glass = PYAFButton(self, "glass", "Define as substrate")
        self.BT_planefit = PYAFButton(self, "planefit", "Planefit")
        self.BT_planefit.setToolTip(
            "Fit a plane (offset and tilt) on the selected substrate ROI and\n"
            "subtract it from the piezo height and topography maps of this\n"
            "file: the substrate becomes the zero height. Only the maps are\n"
            "corrected, not the computations. Define the ROI as substrate\n"
            "first.")
        self.BT_discard = PYAFButton(self, "discard", "Discard data")
        self.BT_discard.setToolTip(
            "Discard data: the curves of the selected ROI are not used in\n"
            "the results (crossed on the maps). The ROI is kept.\n"
            "Keep data: use these curves again.")
        self.BT_volume = PYAFButton(self, "volume", "Compute volume")
        self.BT_volume.setToolTip(
            "Volume between the zero-force topography and the plane of the\n"
            "selected substrate ROI, over the whole map. Define the ROI as\n"
            "substrate first.")

        self.BT_cell = PYAFButton(self, "cell", "Cell geometry")
        self.BT_cell.setToolTip(
            "Plant cells: size, axes and curvatures of the cell drawn as the\n"
            "selected ROI (quadratic surface fitted on its heights), at its\n"
            "centre and at the single curves recorded inside it.")

        # Pencil, below the map with the tools
        self.HL_pencil = QtWidgets.QHBoxLayout()
        self.HL_pencil.addWidget(self.radiobutton_square)
        self.HL_pencil.addWidget(self.radiobutton_circle)
        self.HL_pencil.addWidget(self.radiobutton_dot)
        self.HL_pencil.addWidget(self.input_size)
        self.HL_pencil.addStretch(1)
        self.VL_map.insertLayout(4, self.HL_pencil)

        # ROI assistance, below the pencil (shown with ROI assistance)
        self.CB_auto_mode = PYAFComboBox(self, "auto_mode")
        self.CB_auto_mode.addItem("Regions", "regions")
        self.CB_auto_mode.addItem("Dome", "dome")
        self.CB_auto_mode.fit_to_contents()
        self.CB_auto_mode.setToolTip(
            "Regions: the pixels sorted into the ROIs from all the maps,\n"
            "at least two ROIs marked (cell and glass, domains of a\n"
            "bilayer).\n"
            "Dome: each ROI takes the whole domes of the height its marks\n"
            "touch, separated by the grooves (cells of a tissue); one ROI\n"
            "is enough.")
        self.SL_auto = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.SL_auto.setRange(30, 99)
        self.SL_auto.setValue(70)
        self.SL_auto.setFixedWidth(100)
        self.SL_auto.setToolTip(
            "Pixels sorted with a lower confidence are left uncertain (in\n"
            "no ROI): borders, debris")
        self.label_threshold = QtWidgets.QLabel()
        self.SB_auto_size = QtWidgets.QSpinBox()
        self.SB_auto_size.setRange(1, 9999)
        self.SB_auto_size.setValue(5)
        self.SB_auto_size.setMaximumWidth(70)
        self.SB_auto_size.setToolTip(
            "Groups of fewer pixels of a ROI are left uncertain")
        self.SL_auto.valueChanged.connect(self.auto_classify)
        self.SB_auto_size.valueChanged.connect(self.auto_classify)
        self.SB_dome = QtWidgets.QSpinBox()
        self.SB_dome.setRange(1, 50)
        self.SB_dome.setValue(5)
        self.SB_dome.setSuffix(" %")
        self.SB_dome.setToolTip(
            "Height of a dome above the grooves around it, in % of the\n"
            "range of the heights of the map: lower domes (noise, small\n"
            "bumps) are merged with their neighbors. Lower it if two domes\n"
            "are merged, raise it if a dome is cut")
        self.SB_dome.valueChanged.connect(self.auto_classify)
        self.BT_auto_accept = PYAFButton(self, "auto_accept", "Accept")
        self.set_auto_status("")
        self.BT_auto_clear = PYAFButton(self, "auto_clear", "Clear marks")
        # Options of the regions, or of the domes
        self.W_regions = QtWidgets.QWidget()
        HL_regions = QtWidgets.QHBoxLayout(self.W_regions)
        HL_regions.setContentsMargins(0, 0, 0, 0)
        HL_regions.addWidget(QtWidgets.QLabel("Uncertain below"))
        HL_regions.addWidget(self.SL_auto)
        HL_regions.addWidget(self.label_threshold)
        HL_regions.addSpacing(12)
        HL_regions.addWidget(QtWidgets.QLabel("Min. size [px]"))
        HL_regions.addWidget(self.SB_auto_size)
        self.W_dome = QtWidgets.QWidget()
        HL_dome = QtWidgets.QHBoxLayout(self.W_dome)
        HL_dome.setContentsMargins(0, 0, 0, 0)
        HL_dome.addWidget(QtWidgets.QLabel("Min. dome height"))
        HL_dome.addWidget(self.SB_dome)
        # Two rows, as wide as the map at most (the window does not grow):
        # the mode and the buttons, then the options of the mode
        HL_auto = QtWidgets.QHBoxLayout()
        HL_auto.addWidget(QtWidgets.QLabel("Sort by"))
        HL_auto.addWidget(self.CB_auto_mode)
        HL_auto.addSpacing(12)
        HL_auto.addWidget(self.BT_auto_accept)
        HL_auto.addWidget(self.BT_auto_clear)
        HL_auto.addStretch(1)
        HL_auto_options = QtWidgets.QHBoxLayout()
        HL_auto_options.addWidget(self.W_regions)
        HL_auto_options.addWidget(self.W_dome)
        HL_auto_options.addStretch(1)
        self.W_auto = QtWidgets.QWidget()
        VL_auto = QtWidgets.QVBoxLayout(self.W_auto)
        VL_auto.setContentsMargins(0, 0, 0, 0)
        VL_auto.addLayout(HL_auto)
        VL_auto.addLayout(HL_auto_options)
        self.W_auto.setVisible(False)
        self.VL_map.insertWidget(5, self.W_auto)
        self.update_auto_panel()

        # Table and actions on the ROIs, right of the map
        GL_buttons = QtWidgets.QGridLayout()
        for index, button in enumerate((
                self.button_add, self.button_remove, self.BT_discard,
                self.button_invert, self.BT_glass, self.BT_planefit,
                self.BT_volume, self.BT_cell)):
            GL_buttons.addWidget(button, index // 2, index % 2)
        self.VL_options.addWidget(self.TW, 1)
        self.VL_options.addLayout(GL_buttons)

        self.HL.addLayout(self.VL_options)
        self.VL.addLayout(self.HL)
        self.setLayout(self.VL)

        self.update_widget()

        self.TW.selectRow(0)

    def tableWidget_clicked(self, row, col):
        """Updates the glass button when changing ROI."""
        # Do nothing with these arguments :
        _ = row
        _ = col

        self.update_GUI("BT_glass")

    def button_clicked(self, button):
        """Called when a button is clicked."""
        data = shared.exp.current_data

        if button == "add":
            newid = len(data.roi_list)

            self.logger.debug("Adding ROI %s", newid)

            data.roi_list.append(new_roi(data, newid))
            # The new ROI is the one drawn (Draw column)
            data.roi_selected_row = newid

            self.update_widget()
            self.update_results()

            self.TW.selectRow(data.roi_selected_row)

            reset_results_table()

        elif button == "remove":
            indexes = []
            for item in self.TW.selectedIndexes():
                indexes.append(item.row())

            # It depends which way you have selected so I need to sort the
            # index array
            indexes = numpy.sort(numpy.unique(indexes), -1)[::-1].tolist()

            # Check if there are no glas layers defined in the VTK widget
            # The ROI has to be from the same file and have the same id
            u_id = data.unique_id
            if widgets_list.widget_vtk is not None:
                for layer in shared.layer_list:
                    if layer.type == "glass" and layer.afm_id == u_id:
                        for index in indexes:
                            if index == layer.roi_id:
                                # Display an error message and abort removing
                                ask_user_to_delete_layer_first()
                                self.logger.debug(
                                    "Abort removing ROI %s",
                                    index)
                                return False

            self.logger.debug("Removing ROI %s", indexes)

            if indexes:
                self.update_planefit_roi(data, indexes)
                self.update_substrate_roi(data, indexes)
                for i in indexes:
                    del data.roi_list[i]
                    self.TW.removeRow(i)

                # ROI used to draw: the same one if it is kept, else the
                # first one
                active = data.roi_selected_row
                if not data.roi_list:
                    data.roi_selected_row = None
                elif active is None or active in indexes:
                    data.roi_selected_row = 0
                else:
                    data.roi_selected_row = active - sum(
                        1 for i in indexes if i < active)

                if data.roi_selected_row is not None:
                    self.TW.selectRow(0)

                self.update_widget()
                if self.auto_active():
                    # Without the marks of the ROIs removed
                    self.auto_update()
                self.rois_changed(data)
                widgets_list.widget_compute.update_GUI("roi_list_height_corr")

            reset_results_table()

        elif button == "bg":
            value = self.BG.checkedId()
            if value == 0:
                shape = "square"
            elif value == 1:
                shape = "circle"
            elif value == 2:
                shape = "dot"
            shared.exp.roi_cursor_shape = shape

            self.logger.debug("Changing roi shape to %s", shape)

        elif button == "tool":
            self.change_cursor_mode(
                ("select", "roi_add", "roi_erase")[self.BG_tool.checkedId()])

        elif button == "invert":
            for item in self.TW.selectionModel().selectedRows():
                index = item.row()

                newid = len(data.roi_list)

                txt = (
                    "Adding ROI (by inverting). New ROI id = %s, "
                    "old ROI = %s")
                self.logger.debug(txt, newid, index)

                data.roi_list.append(new_roi(data, newid))

                # Get new values (can be slow)
                inverted_roi = []
                for i in range(data.nbr_pixels_x):
                    for j in range(data.nbr_pixels_y):
                        values = data.roi_list[index].values
                        if not math_tools.in_list(values, [i, j]):
                            inverted_roi.append([i, j])

                # Add data to ROI
                data.roi_list[newid].values = inverted_roi

            self.update_widget()
            self.rois_changed(data)
            reset_results_table()

        elif button == "button_group_selected":
            data.roi_selected_row = self.BG_selected.checkedId()
            self.logger.debug("Selecting ROI %s", data.roi_selected_row)

        elif button == "glass":
            for item in self.TW.selectionModel().selectedRows():
                index = item.row()
                roi = data.roi_list[index]
                if roi.glass_coeffs is None:
                    self.logger.debug("Set as glass. ROI %s", index)
                    # Plane of the substrate (fitted again on the
                    # topography by the elasticity)
                    roi.glass_coeffs = planefit.fit_substrate(data, roi)
                else:
                    self.logger.debug("Remove glass. ROI %s", index)

                    roi.glass_coeffs = None

                    # The planefit used this substrate
                    if data.planefit_applied and data.planefit_roi == index:
                        self.set_planefit(data, None)
                bec.mark_outdated(data)

            self.update_widget()
            widgets_list.widget_compute.update_GUI("roi_list_height_corr")
            widgets_list.widget_compute.update_GUI("BT_height_correction")
            # The substrate is not a cell
            plant.refresh_results()

        elif button == "volume":
            self.show_volume(data)

        elif button == "planefit":
            rows = [item.row()
                    for item in self.TW.selectionModel().selectedRows()]
            if planefit.get_planefit_roi(data) is not None:
                # Remove the planefit
                self.set_planefit(data, None)
            elif rows:
                self.set_planefit(data, rows[0])
            self.update_GUI("BT_glass")

        elif button == "cell":
            self.show_cell_geometry(data)

        elif button == "auto_accept":
            self.auto_accept(data)

        elif button == "auto_clear":
            self.auto_marks.pop(data.unique_id, None)
            self.auto_erased.pop(data.unique_id, None)
            self.auto_update()
            refresh_rois()

        elif button == "discard":
            # Discard (or keep again) the curves of the selected ROI
            rows = [item.row()
                    for item in self.TW.selectionModel().selectedRows()]
            if not rows:
                return
            roi = data.roi_list[rows[0]]
            if getattr(roi, "discarded", []):
                # Keep data: only the curves discarded by this ROI (not the
                # ones discarded otherwise)
                for i, j in roi.discarded:
                    data.discarded_curves[i][j] = 0
                roi.discarded = []
            else:
                roi.discarded = []
                for i, j in roi.values:
                    if data.discarded_curves[i][j] == 0:
                        data.discarded_curves[i][j] = 1
                        roi.discarded.append([i, j])
            self.update_widget()

            # Refresh the meshgrids and the results
            widgets_list.widget_data.meshgrid.update_plot()
            widgets_list.widget_results.update_MPL("MPL_canvas1")
            self.update_map()
            self.update_results()

    def show_cell_geometry(self, data):
        """Geometry of the cell drawn as the selected ROI (see
        plant_tools.cell_geometry), at its centre and at the single curves
        recorded inside it."""
        rows = [item.row() for item in self.TW.selectionModel().selectedRows()]
        if not rows:
            return
        roi = data.roi_list[rows[0]]
        inside = {(int(i), int(j)) for i, j in roi.values}
        curves = [(index, position) for index, _, position
                  in plant_tools.curves_on_map(data)
                  if (position["i"], position["j"]) in inside]
        geometry = plant_tools.cell_geometry(
            data, roi.values, [[p["x"], p["y"]] for _, p in curves])
        title = "Cell geometry - ROI " + str(rows[0] + 1)
        if geometry is None:
            QtWidgets.QMessageBox.warning(
                self, title, "Draw the cell as a ROI first (at least 6 "
                "pixels).")
            return

        def number(value, digits=4):
            return "%.*g" % (digits, value)

        names = ["Centre"] + ["File " + str(index + 1) for index, _ in curves]
        rows_html = ""
        saddle = False
        for name, point in zip(names, geometry["points"]):
            saddle |= not point["dome"]
            radii = ["%.1f" % (1.0 / value) if value > 0 else "-"
                     for value in (point["k_major"], point["k_minor"])]
            rows_html += (
                "<tr><td>" + name + "</td><td>" +
                number(point["x"] / 1000.0, 3) + ", " +
                number(point["y"] / 1000.0, 3) + "</td><td>" +
                number(point["k_major"] * 1e3) + "</td><td>" +
                number(point["k_minor"] * 1e3) + "</td><td>" +
                radii[0] + ", " + radii[1] + "</td><td>" +
                number(point["kM"] * 1e3) + "</td><td>" +
                number(point["kG"] * 1e6) + "</td><td>" +
                number(point["f"], 3) + "</td><td>" +
                ("dome" if point["dome"] else "<b>saddle</b>") +
                "</td></tr>")
        text = (
            "<p>" + data.filename + "<br>Heights: " + geometry["heights"] +
            " (without flatten)</p>"
            "<p>Area: " + number(geometry["area"]) + " µm² (" +
            str(geometry["pixels"]) + " pixels)<br>"
            "Major axis: " + number(geometry["angle"], 3) + "° from x<br>"
            "Half-lengths (centre to edge): L<sub>major</sub> = " +
            number(geometry["L_major"]) + " µm, L<sub>minor</sub> = " +
            number(geometry["L_minor"]) + " µm<br>"
            "Quadratic surface fit: RMS " + number(geometry["rms"], 3) +
            " nm</p>"
            "<table cellspacing='6'><tr><th></th><th>x, y [µm]</th>"
            "<th>κ<sub>major</sub> [1/mm]</th><th>κ<sub>minor</sub> [1/mm]"
            "</th><th>Radii [µm]</th><th>κ<sub>M</sub> [1/mm]</th>"
            "<th>κ<sub>G</sub> [1/mm²]</th><th>f</th><th></th></tr>" +
            rows_html + "</table>")
        if saddle:
            text += ("<p>Saddle: the surface is not a dome there (junction "
                     "between cells?), the shell models do not apply.</p>")
        box = QtWidgets.QMessageBox(self)
        box.setWindowTitle(title)
        box.setText(title)
        box.setInformativeText(text)
        box.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        box.exec()

    def set_planefit(self, data, roi_index):
        """Apply the planefit with a substrate ROI (None: remove it), then
        rescale and redraw the height maps."""
        self.logger.debug("Planefit with ROI %s", roi_index)
        data.planefit_applied = roi_index is not None
        data.planefit_roi = roi_index
        heights_changed([data])

    def show_volume(self, data):
        """Computes the volume above the selected substrate ROI and displays
        it."""
        rows = [item.row() for item in self.TW.selectionModel().selectedRows()]
        if not rows or data.roi_list[rows[0]].glass_coeffs is None:
            return
        res = planefit.get_volume(data, data.roi_list[rows[0]])
        x_size, y_size = res["pixel_size"]
        title = ("Volume relative to the plane of the substrate ROI " +
                 str(rows[0] + 1))
        text = ("Zero-force topography, whole map\n\n"
                "Net volume: %.4g µm³\n"
                "Above the plane: %.4g µm³\n"
                "Below the plane: %.4g µm³\n"
                "Maximum height: %.4g nm\n\n"
                "%d pixels of %.4g × %.4g nm") % (
                    res["net"] * 1e-9, res["above"] * 1e-9,
                    res["below"] * 1e-9, res["max_height"], res["pixels"],
                    x_size, y_size)
        if res["no_topography"]:
            text += ("\n%d pixels without topography (no contact point) are "
                     "not counted") % res["no_topography"]
        box = QtWidgets.QMessageBox(self)
        box.setWindowTitle("Volume - " + data.filename)
        box.setText(title)
        box.setInformativeText(text)
        box.setTextInteractionFlags(QtCore.Qt.TextSelectableByMouse)
        box.exec()

    @staticmethod
    def update_3D(data):
        """Redraw the 3D layers of the file (tomography and substrates)."""
        if widgets_list.widget_vtk is None:
            return
        for layer in shared.layer_list:
            if layer.type not in ("afm", "glass"):
                continue
            if shared.exp.list[layer.afm_id] is not data:
                continue
            if layer.type == "afm":
                layer.update_tomography()
            else:
                layer.update_actor()

    def update_planefit_roi(self, data, removed):
        """Keep the index of the planefit ROI when ROIs are removed (no
        planefit if its ROI is removed)."""
        if data.planefit_roi is None:
            return
        if data.planefit_roi in removed:
            data.planefit_applied = False
            data.planefit_roi = None
            planefit.reset_color_scales(data)
            self.update_3D(data)
        else:
            data.planefit_roi -= sum(1 for i in removed
                                     if i < data.planefit_roi)

    @staticmethod
    def update_substrate_roi(data, removed):
        """Keep the substrate ROI of the bottom effect correction when ROIs
        are removed (roi_glass_id: its position + 1; none if it is
        removed)."""
        if not data.roi_glass_id:
            return
        index = data.roi_glass_id - 1
        if index in removed:
            data.roi_glass_id = 0
        else:
            data.roi_glass_id -= sum(1 for i in removed if i < index)

    def update_results(self):
        """Update the results in the results tables."""
        main = widgets_list.widget_main
        if main.tabs.currentIndex() in (3, 4):
            # A plots tab is displayed: recompute the values now (the values
            # are recomputed when going to the plots tabs from the others)
            widgets_list.widget_results_single.schedule_refresh()
        else:
            main.update_results = True

    def rois_changed(self, data, pixels=True):
        """The ROIs of a file changed (pixels, or display and colors with
        pixels=False): maps redrawn, and with new pixels the substrate
        planes fitted again and the values which depend on the ROIs (B.E.C.
        out of date, plant cells, results tables)."""
        refresh_rois()
        if not pixels:
            return
        substrate = data.roi_list[data.roi_glass_id - 1] \
            if 0 < data.roi_glass_id <= len(data.roi_list) else None
        before = None if substrate is None or \
            substrate.glass_coeffs is None else \
            numpy.array(substrate.glass_coeffs, dtype=float)
        planefit.update_substrates(data)
        if before is not None and not numpy.allclose(
                before, numpy.array(substrate.glass_coeffs, dtype=float)):
            # New plane of the substrate of the B.E.C.
            bec.mark_outdated(data)
            widgets_list.widget_compute.update_GUI("BT_height_correction")
        if planefit.get_planefit_roi(data) is not None:
            planefit.reset_color_scales(data)
        plant.refresh_results()
        self.update_results()

    def checkbox_clicked(self, checkbox, cb_id):
        """Called when a checkbox is clicked."""
        data = shared.exp.current_data

        if checkbox == "display_roi":
            value = self.TW.cellWidget(cb_id, 2).checkbox.isChecked()
            self.logger.debug("Display/hide ROI %s, %s", cb_id, value)
            data.roi_list[cb_id].display = value
            self.rois_changed(data, pixels=False)

    def update_lists(self):
        """Lists of the files and of the maps computed for the current
        file; the map is drawn again for another file."""
        data = shared.exp.current_data
        self.list_files.blockSignals(True)
        self.list_files.clear()
        for index, dt in enumerate(shared.exp.list):
            if not dt.is_single:
                self.list_files.addItem(str(index + 1) + " : " +
                                        dt.filename, index)
        position = self.list_files.findData(shared.exp.id_selected)
        if position >= 0:
            self.list_files.setCurrentIndex(position)
        self.list_files.blockSignals(False)

        self.list_maps.blockSignals(True)
        self.list_maps.clear()
        maps = map_values.available_maps(data) if not data.is_single else []
        for map_type in maps:
            self.list_maps.addItem(map_values.MAPS[map_type][0], map_type)
        shown = shared.exp.roi_meshgrid_type
        self.list_maps.setCurrentIndex(max(self.list_maps.findData(
            shown if shown in maps else "piezo"), 0))
        self.list_maps.fit_to_contents()
        self.list_maps.blockSignals(False)

        # Single curve: no map, no ROI
        for widget in (self.W_canvas, self.RB_select, self.RB_draw,
                       self.RB_erase, self.list_maps):
            widget.setEnabled(not data.is_single)
        self.TG_auto.setEnabled(not data.is_single and auto_roi.available())
        self.update_auto_panel()
        self.label_no_map.setVisible(data.is_single)
        if data.is_single and self.cursor_mode != "select":
            # Tool of the maps, back with the next map
            self.map_tool = self.cursor_mode
            self.RB_select.setChecked(True)
            self.change_cursor_mode("select")
        elif not data.is_single and self.map_tool is not None:
            self.BG_tool.button(("select", "roi_add", "roi_erase").index(
                self.map_tool)).setChecked(True)
            self.change_cursor_mode(self.map_tool)
            self.map_tool = None

        if self.shown_file is not data:
            self.shown_file = data
            # Proposal of this file (from its marks)
            self.auto_features = None
            self.auto_update()
            self.update_map()

    def update_map(self):
        """Draw the map of the window again (empty for a single curve)."""
        # The maps may have changed (computations, planefit, flatten): the
        # features of ROI assistance are computed again at the next proposal
        self.auto_features = None
        self.auto_domes = None
        self.MPL_canvas.update_plot()
        if shared.exp.current_data.is_single:
            # No map: no plot to draw the ROIs on
            self.MPL_canvas.canvas = None
            self.MPL_canvas.draw_idle()

    def change_cursor_mode(self, mode):
        """Tool used on the map: select a curve, draw or erase a ROI."""
        self.cursor_mode = mode
        if mode == "select":
            self.MPL_canvas.unsetCursor()
        else:
            self.MPL_canvas.setCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
            if isinstance(self.MPL_canvas.canvas, MainPlot):
                self.MPL_canvas.canvas.first_blit = True

    def map_plot(self, event):
        """Plot of the map if the event is on it, else None."""
        plot = self.MPL_canvas.canvas
        if isinstance(plot, MainPlot) and event.inaxes == plot.axes:
            return plot
        return None

    def canvas_press(self, event):
        """Click on the map: select a curve, or start drawing."""
        data = shared.exp.current_data
        if event.button != 1 or data.is_single or \
                self.map_plot(event) is None:
            return
        xpos, ypos = misc_tools.get_position_on_meshgrid(event)
        if self.cursor_mode == "select":
            widgets_list.widget_main.change_curve(xpos, ypos)
            refresh_rois()
        elif data.roi_selected_row is None:
            QtWidgets.QMessageBox.information(
                self, "ROI Manager", "Add a ROI first (Add ROI).")
        else:
            self.drawing = True
            self.draw_roi(xpos, ypos)

    def canvas_motion(self, event):
        """Mouse over the map: draw (button pressed) or preview the
        pencil."""
        data = shared.exp.current_data
        if self.cursor_mode == "select" or data.is_single or \
                data.roi_selected_row is None or \
                self.map_plot(event) is None:
            return
        xpos, ypos = misc_tools.get_position_on_meshgrid(event)
        if self.drawing and event.button == 1:
            self.draw_roi(xpos, ypos)
        else:
            self.pre_draw_roi(xpos, ypos)

    def canvas_release(self, event):
        """End of the drawing (also outside of the map)."""
        _ = event
        if not self.drawing:
            return
        self.drawing = False
        data = shared.exp.current_data
        if self.auto_active():
            # Marks changed: new proposal (the ROIs do not change)
            self.auto_update()
            refresh_rois()
            return
        self.update_widget()
        self.rois_changed(data)
        reset_results_table()

    def pre_draw_roi(self, xpos, ypos):
        """Pencil drawn on the map, without changing the ROI."""
        data = shared.exp.current_data
        local_roi, local_roi_mask = prepare_roi(xpos, ypos)
        local_roi, local_roi_mask = \
            check_roi_boundaries(data, local_roi, local_roi_mask)
        mode = "preview_add" if self.cursor_mode == "roi_add" \
            else "preview_erase"
        self.MPL_canvas.canvas.update_blit(mode, local_roi, local_roi_mask)

    def draw_roi(self, xpos, ypos):
        """Pixels under the pencil added to (or removed from) the ROI
        checked in the Draw column."""
        data = shared.exp.current_data
        roi = data.roi_list[data.roi_selected_row]
        local_roi, local_roi_mask = prepare_roi(xpos, ypos)
        local_roi, local_roi_mask = \
            check_roi_boundaries(data, local_roi, local_roi_mask)
        if self.auto_active():
            # Marks of the ROI (with ROI assistance), not its pixels: Draw
            # marks pixels in the ROI, Erase takes pixels out of it (by
            # dome, the region flooded from them, see auto_classify)
            marks = self.auto_marks.setdefault(data.unique_id, {}) \
                .setdefault(roi, set())
            erased = self.auto_erased.setdefault(data.unique_id, {}) \
                .setdefault(roi, set())
            pencil = {tuple(pixel) for index, pixel in enumerate(local_roi)
                      if local_roi_mask[index] and pixel is not None}
            if self.cursor_mode == "roi_add":
                marks |= pencil
                erased -= pencil
                mode = "local_add"
            else:
                marks -= pencil
                erased |= pencil
                mode = "local_erase"
            self.MPL_canvas.canvas.update_blit(mode, local_roi,
                                               local_roi_mask)
            return
        inside = {tuple(value) for value in roi.values}
        if self.cursor_mode == "roi_add":
            # Only the new pixels are drawn (transparent color)
            added = []
            for index, pixel in enumerate(local_roi):
                new = local_roi_mask[index] and pixel is not None and \
                    tuple(pixel) not in inside
                if new:
                    roi.values.append(pixel)
                    inside.add(tuple(pixel))
                added.append(new)
            local_roi_mask = added
            mode = "local_add"
        else:
            erased = {tuple(pixel) for index, pixel in enumerate(local_roi)
                      if local_roi_mask[index] and pixel is not None}
            roi.values[:] = [value for value in roi.values
                             if tuple(value) not in erased]
            mode = "local_erase"
        self.MPL_canvas.canvas.update_blit(mode, local_roi, local_roi_mask)

    # --- ROI assistance ------------------------------------------------------

    def auto_active(self):
        """True if the tools mark the pixels of the ROIs (ROI assistance)."""
        return self.TG_auto.isChecked() and \
            not shared.exp.current_data.is_single

    def toggle_clicked(self, name):
        """Called when a toggle is clicked."""
        if name == "roi_assistance":
            self.auto_toggled(self.TG_auto.isChecked())

    def auto_toggled(self, checked):
        """ROI assistance on: proposal from the marks of the file (the pixels
        of the ROIs already drawn for the ROIs without marks); off: the ROIs
        drawn again (the marks are kept for the session)."""
        data = shared.exp.current_data
        if checked and not data.is_single:
            # The ROIs already drawn are the marks of the ROIs without marks
            marks = self.auto_marks.setdefault(data.unique_id, {})
            for roi in data.roi_list:
                if roi.values and not marks.get(roi):
                    marks[roi] = {(int(i), int(j)) for i, j in roi.values}
        self.update_auto_panel()
        self.auto_update()
        refresh_rois()

    def update_auto_panel(self):
        """Options of ROI assistance, shown with it."""
        self.W_auto.setVisible(self.auto_active())
        by_dome = self.auto_by_dome()
        # The options of the other mode hidden first: both in the row would
        # enlarge the window (not shrunk back)
        hidden, shown = (self.W_regions, self.W_dome) if by_dome else \
            (self.W_dome, self.W_regions)
        hidden.setVisible(False)
        shown.setVisible(True)
        self.label_threshold.setText("%d %%" % self.SL_auto.value())
        result = self.auto_result
        self.BT_auto_accept.setEnabled(
            result is not None and result["classes"] is not None)
        data = shared.exp.current_data
        self.BT_auto_clear.setEnabled(any(
            self.auto_marks.get(data.unique_id, {}).values()) or any(
            self.auto_erased.get(data.unique_id, {}).values()))

    def auto_by_dome(self):
        """True if the pixels are sorted by dome (domes of the height),
        False by regions (forest)."""
        return self.CB_auto_mode.currentData() == "dome"

    def auto_mode_changed(self):
        """Sorted by regions or by dome: new proposal."""
        self.auto_update()
        refresh_rois()

    def auto_dome_labels(self, data, seeds):
        """Domes of the file, and the regions flooded from seeds (see
        auto_roi.domes), computed again when they change."""
        key = (data.unique_id, self.SB_dome.value(),
               tuple(frozenset(pixels) for pixels in seeds))
        if self.auto_domes is None or self.auto_domes[0] != key:
            self.auto_domes = (key, auto_roi.domes(
                auto_roi.height(data), self.SB_dome.value() / 100.0, seeds))
        return self.auto_domes[1]

    def set_auto_status(self, text):
        """State of the proposal (pixels of each ROI, or what is missing),
        in the tooltip of Accept."""
        self.auto_status = text
        tip = "The pixels of the proposal become the pixels of the ROIs\n" \
            "marked, and ROI assistance is turned off (the ROIs can be\n" \
            "edited with Draw and Erase)."
        if text:
            tip += "\n\n" + text
        self.BT_auto_accept.setToolTip(tip)

    def auto_rois(self, data):
        """ROIs of the file with marks, in the order of the list, and their
        marks."""
        marks = self.auto_marks.get(data.unique_id, {})
        erased = self.auto_erased.get(data.unique_id, {})
        # Not the ROIs removed
        for store in (marks, erased):
            for roi in [roi for roi in store if roi not in data.roi_list]:
                del store[roi]
        rois = [roi for roi in data.roi_list if marks.get(roi)]
        return rois, [marks[roi] for roi in rois]

    def auto_update(self):
        """Proposal from the marks of the file of the map: by regions, the
        forest learns the ROIs on the marked pixels and gives the
        probabilities of the ROIs for all the pixels; by dome, the domes
        under the marks (see auto_classify and auto_roi)."""
        self.auto_result = None
        if not self.auto_active():
            self.update_auto_panel()
            return
        data = shared.exp.current_data
        rois, marks = self.auto_rois(data)
        shape = (data.nbr_pixels_x, data.nbr_pixels_y)
        labels = auto_roi.marks_labels(shape, marks)
        self.auto_result = {"file": data.unique_id, "rois": rois,
                            "marks": labels, "proba": None, "classes": None}
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.WaitCursor)
        try:
            if self.auto_by_dome():
                if not rois:
                    raise auto_roi.AutoRoiError(
                        "Mark the domes (a stroke or a click in each) with "
                        "the ROIs (Draw column).")
            elif len(rois) < 2:
                raise auto_roi.AutoRoiError(
                    "Mark pixels in at least two ROIs (Draw column) to get "
                    "a proposal.")
            if self.auto_features is None or \
                    self.auto_features[0] != data.unique_id:
                self.auto_features = (data.unique_id,
                                      auto_roi.features(data))
            self.auto_result["proba"] = auto_roi.probabilities(
                self.auto_features[1], labels, data.inside_region())
        except auto_roi.AutoRoiError as error:
            self.set_auto_status(str(error))
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.auto_classify()

    def auto_classify(self):
        """ROIs of the pixels: by regions from the probabilities, with the
        confidence threshold and the minimum size (fast: not learned
        again); by dome, the domes under the marks of each ROI."""
        result = self.auto_result
        data = shared.exp.current_data
        if result is None or not result["rois"]:
            self.update_auto_panel()
            return
        erased = self.auto_erased.get(data.unique_id, {})
        if self.auto_by_dome():
            # The pixels erased from a ROI are flooded as regions it does
            # not take (e.g. a flat substrate next to a cell, which has no
            # dome of its own)
            seeds = [erased.get(roi, set()) for roi in result["rois"]]
            try:
                dome_labels, number = self.auto_dome_labels(data, seeds)
            except auto_roi.AutoRoiError as error:
                self.set_auto_status(str(error))
                self.update_auto_panel()
                return
            result["classes"], several = auto_roi.dome_classes(
                dome_labels, result["marks"], data.inside_region(),
                {index: {number + 1 + index}
                 for index in range(len(seeds))})
            end = "%d dome%s in the map" % (number, "s" if number > 1
                                            else "")
            if several:
                end += " (%d marked in several ROIs: left out; with no " \
                    "dome, sort by Regions)" % several
        elif result["proba"] is not None:
            result["classes"] = auto_roi.classify(
                result["proba"], result["marks"],
                self.SL_auto.value() / 100.0, self.SB_auto_size.value())
            end = None
        else:
            self.update_auto_panel()
            return
        # Pixels erased from the ROIs
        for index, roi in enumerate(result["rois"]):
            for i, j in erased.get(roi, ()):
                if result["classes"][i, j] == index:
                    result["classes"][i, j] = -1
        if end is None:
            end = "uncertain: %d px" % numpy.count_nonzero(
                result["classes"] == -1)
        counts = ["ROI %d: %d px" % (data.roi_list.index(roi) + 1,
                                     numpy.count_nonzero(
                                         result["classes"] == index))
                  for index, roi in enumerate(result["rois"])]
        self.set_auto_status(", ".join(counts + [end]))
        if self.sender() in (self.SL_auto, self.SB_auto_size, self.SB_dome):
            refresh_rois()
        self.update_auto_panel()

    def auto_overlay(self, data):
        """What the map of the window shows with ROI assistance, None without:
        the marks of the ROIs ([pixels, color]), and the proposal (ROI of
        each pixel, see auto_roi.classify, and the colors of the ROIs) or
        None."""
        if not self.auto_active() or data is not shared.exp.current_data:
            return None
        rois, marks = self.auto_rois(data)
        colors = [data.roi_color_names[roi.color] for roi in rois]
        result = self.auto_result
        classes = None
        if result is not None and result["file"] == data.unique_id and \
                result["classes"] is not None:
            classes = result["classes"]
            colors = [data.roi_color_names[roi.color]
                      for roi in result["rois"]]
        return {"marks": [([list(pixel) for pixel in pixels],
                           data.roi_color_names[roi.color])
                          for roi, pixels in zip(rois, marks)],
                "classes": classes, "colors": colors}

    def auto_accept(self, data):
        """The pixels of the proposal become the pixels of the ROIs marked
        (after a confirmation if pixels of the ROIs are not in it)."""
        result = self.auto_result
        if result is None or result["classes"] is None or \
                result["file"] != data.unique_id:
            return
        rois = [roi for roi in result["rois"] if roi in data.roi_list]
        classes = result["classes"]
        reduced = [str(data.roi_list.index(roi) + 1)
                   for index, roi in enumerate(result["rois"])
                   if roi in data.roi_list and
                   any(classes[i, j] != index for i, j in roi.values)]
        if reduced:
            answer = QtWidgets.QMessageBox.question(
                self, "ROI Manager",
                "Pixels of the ROI%s %s are not in the proposal. Replace "
                "the pixels with the proposal?" %
                ("s" if len(reduced) > 1 else "", ", ".join(reduced)))
            if answer != QtWidgets.QMessageBox.Yes:
                return
        for index, roi in enumerate(result["rois"]):
            if roi in data.roi_list:
                roi.values = [[int(i), int(j)] for i, j in
                              numpy.argwhere(result["classes"] == index)]
        self.logger.debug("ROI assistance accepted for ROIs %s", [
            data.roi_list.index(roi) + 1 for roi in rois])
        # The ROIs, editable with the tools
        self.TG_auto.setChecked(False)
        self.auto_toggled(False)
        self.update_widget()
        self.rois_changed(data)
        reset_results_table()

    def update_widget(self):
        """Update the content of the ROI manager."""
        self.update_lists()

        # Get previously selected rows
        selected_rows = []
        for item in self.TW.selectionModel().selectedRows():
            selected_rows.append(item.row())

        data = shared.exp.current_data
        for i in range(self.TW.rowCount(), -1, -1):
            self.TW.removeRow(i)

        self.BG_selected = PYAFButtonGroup(self, "button_group_selected")

        # Add the rows to the tablewidget
        for i in range(len(data.roi_list)):
            self.TW.insertRow(i)
            self.TW.setRowHeight(i, 32)

            # Index
            label = QtWidgets.QLabel()
            label.setFont(widgets_list.widget_main.smallfont)
            label.setText(str(i + 1))
            self.TW.setCellWidget(i, 0, gui_tools.CenteredCellwidget(label))

            # Radiobutton
            self.TW.setCellWidget(i, 1, gui_tools.CenteredCellRadioButton())
            button = self.TW.cellWidget(i, 1).radio_button
            self.BG_selected.addButton(button, i)
            if data.roi_selected_row == i or data.roi_selected_row is None:
                button.setChecked(True)
                data.roi_selected_row = i

            # Checkbox
            cb = gui_tools.CenteredCellCheckbox(self, "display_roi", i)
            self.TW.setCellWidget(i, 2, cb)
            disp = data.roi_list[i].display
            self.TW.cellWidget(i, 2).checkbox.setChecked(disp)

            # Pixels
            label = QtWidgets.QLabel()
            label.setFont(widgets_list.widget_main.smallfont)
            label.setText(str(len(data.roi_list[i].values)))
            self.TW.setCellWidget(i, 3, gui_tools.CenteredCellwidget(label))

            # Area
            xsize = data.scan_size_x / data.nbr_pixels_x
            ysize = data.scan_size_y / data.nbr_pixels_y
            one_area = xsize * ysize
            label = QtWidgets.QLabel()
            label.setFont(widgets_list.widget_main.smallfont)
            area = len(data.roi_list[i].values) * one_area * 1e-6
            label.setText("%.4g" % area)
            self.TW.setCellWidget(i, 4, gui_tools.CenteredCellwidget(label))

            # Color
            colorlist = PYAFComboBox(self, "colorlist")
            colorlist.setFont(widgets_list.widget_main.smallfont)
            for j in range(len(data.roi_color_names)):
                colorlist.addItem(data.roi_color_names[j])
            colorlist.setCurrentIndex(data.roi_list[i].color)
            # In a layout, to keep the size of the cell
            self.TW.setCellWidget(i, 5, gui_tools.FullCellWidget(colorlist))

            # Glass
            label = QtWidgets.QLabel()
            label.setFont(widgets_list.widget_main.smallfont)
            if data.roi_list[i].glass_coeffs is not None:
                label.setText("Yes")
            else:
                label.setText("")
            self.TW.setCellWidget(i, 6, gui_tools.CenteredCellwidget(label))

            # Discarded data
            label = QtWidgets.QLabel()
            label.setFont(widgets_list.widget_main.smallfont)
            if getattr(data.roi_list[i], "discarded", []):
                label.setText("Discarded")
            self.TW.setCellWidget(i, 7, gui_tools.CenteredCellwidget(label))

        # Select rows
        for row in selected_rows:
            self.TW.selectRow(row)

        # Select the first row
        if selected_rows == []:
            self.TW.selectRow(0)

        # Shape
        if shared.exp.roi_cursor_shape == "square":
            self.radiobutton_square.setChecked(True)
        elif shared.exp.roi_cursor_shape == "circle":
            self.radiobutton_circle.setChecked(True)
        elif shared.exp.roi_cursor_shape == "dot":
            self.radiobutton_dot.setChecked(True)

        # Size
        self.input_size.changeValue(shared.exp.roi_cursor_size)

        self.update_GUI("BT_glass")

    def update_GUI(self, element):
        """Update some elements of the ROI manager."""
        data = shared.exp.current_data

        if element == "BT_glass":
            found_item = False
            for item in self.TW.selectionModel().selectedRows():
                index = item.row()
                roi = data.roi_list[index]
                found_item = True

                # Before the elasticity: on the piezo height (fitted again
                # on the topography by the elasticity)
                self.BT_glass.setEnabled(bool(roi.values) and
                                         not data.is_single)

                if roi.glass_coeffs is None:
                    self.BT_glass.setText("Define as substrate")
                else:
                    self.BT_glass.setText("No substrate")

        if not found_item:
            # No Roi defined, disable glass button
            self.BT_glass.setEnabled(False)

        if element == "BT_glass":
            # Planefit: with the selected ROI if it is a substrate, or to
            # remove the current planefit
            applied = planefit.get_planefit_roi(data) is not None
            rows = [item.row()
                    for item in self.TW.selectionModel().selectedRows()]
            substrate = bool(rows) and \
                data.roi_list[rows[0]].glass_coeffs is not None
            self.BT_planefit.setText(
                "Remove planefit" if applied else "Planefit")
            self.BT_planefit.setEnabled(applied or substrate)
            self.BT_volume.setEnabled(substrate and data.stiffness_calculated)
            discarded = bool(rows) and \
                bool(getattr(data.roi_list[rows[0]], "discarded", []))
            self.BT_discard.setText("Keep data" if discarded
                                    else "Discard data")
            self.BT_discard.setEnabled(discarded or (
                bool(rows) and bool(data.roi_list[rows[0]].values)))

    def input_updated(self, name):
        """Called when an input is updated."""
        if name == "size":
            value = self.input_size.get_int_value()
            shared.exp.roi_cursor_size = value
            self.logger.debug("Change roi pencil size to %s", value)

    def list_updated(self, name):
        """Called when a list is updated."""
        data = shared.exp.current_data

        if name == "colorlist":
            for i in range(len(data.roi_list)):
                cell = self.TW.cellWidget(i, 5)
                data.roi_list[i].color = cell.widget.currentIndex()
            self.rois_changed(data, pixels=False)

        elif name == "files":
            index = self.list_files.currentData()
            if index is not None and index != shared.exp.id_selected:
                # Same file in the main tabs
                widgets_list.widget_results.list_exp.setCurrentIndex(index)
                widgets_list.widget_main.file_changed("list1")

        elif name == "auto_mode":
            self.auto_mode_changed()

        elif name == "maps":
            map_type = self.list_maps.currentData()
            if map_type is not None:
                shared.exp.roi_meshgrid_type = map_type
                self.update_map()


def new_roi(data, roi_id):
    """New ROI, with the first color not used by the other ROIs of the file
    (the colors are used again when all of them are used)."""
    roi = experiment.ROI(roi_id)
    used = [other.color for other in data.roi_list]
    nbr_colors = len(data.roi_color_names)
    free = [i for i in range(nbr_colors) if i not in used]
    roi.color = free[0] if free else len(data.roi_list) % nbr_colors
    return roi


def heights_changed(datas):
    """The heights of files changed (flatten or planefit): zero-force
    topography with the flatten, substrate planes fitted again (B.E.C. out of
    date), color scales, and the maps, profiles, 3D view and values which
    show the heights."""
    for data in datas:
        if data.sync_topography() and planefit.update_substrates(data) and \
                data.roi_glass_id:
            bec.mark_outdated(data)
        planefit.reset_color_scales(data)
    compute = widgets_list.widget_compute
    if compute is not None:
        compute.update_MPL("MPL_meshgrid")
        compute.update_GUI("BT_height_correction")
    results = widgets_list.widget_results
    if results is not None:
        results.update_MPL("MPL_canvas1")
    if widgets_list.widget_profiles is not None:
        widgets_list.widget_profiles.update_widget()
    manager = widgets_list.widget_roi_manager
    if manager is not None:
        manager.update_map()
    if widgets_list.widget_multimeshgrids is not None:
        widgets_list.widget_multimeshgrids.update_widget()
    for data in datas:
        RoiManagerWidget.update_3D(data)
    # Values of the heights and of the plant cells in the results
    plant.refresh_results()
    main = widgets_list.widget_main
    if main is not None:
        if main.tabs.currentIndex() in (3, 4):
            widgets_list.widget_results_single.schedule_refresh()
        else:
            main.update_results = True


def open_roi_manager():
    """Open the ROI manager (or bring it to the front)."""
    if widgets_list.widget_roi_manager is None:
        RoiManagerWidget(widgets_list.widget_results)
        widgets_list.widget_roi_manager.show()
    else:
        widgets_list.widget_roi_manager.update_widget()
        widgets_list.widget_roi_manager.activateWindow()
        widgets_list.widget_roi_manager.raise_()


def refresh_rois():
    """ROIs drawn again on all the maps showing them (ROI manager, data and
    results tabs, multiple meshgrids)."""
    for plot in (getattr(widgets_list.widget_roi_manager, "MPL_canvas", None),
                 getattr(widgets_list.widget_results, "MPL_canvas1", None),
                 getattr(widgets_list.widget_data, "meshgrid", None)):
        if plot is not None and isinstance(getattr(plot, "canvas", None),
                                           MainPlot):
            plot.canvas.update_blit("all")
    update_roi_on_multimeshgrids_widget()


def reset_results_table():
    """Rows of the ROIs in the results table."""
    if widgets_list.widget_results_single is not None:
        widgets_list.widget_results_single.reset_table()


def prepare_roi(xpos, ypos):
    """Pixels under the pencil at a position of the map (from 1): list of
    [i, j], and mask of the ones in the pencil."""
    shape = shared.exp.roi_cursor_shape
    dot_size = shared.exp.roi_cursor_size
    local_roi = []
    local_roi_mask = []
    if shape == "square":
        for i in range(-dot_size - 1, dot_size):
            for j in range(-dot_size - 1, dot_size):
                local_roi.append([xpos + i, ypos + j])
                local_roi_mask.append(True)
    elif shape == "dot":
        local_roi.append([xpos - 1, ypos - 1])
        local_roi_mask.append(True)
    elif shape == "circle":
        circle_limit = math.sqrt((dot_size - 0.5) ** 2 + (dot_size - 0.5) ** 2)
        for i in range(-dot_size, dot_size + 1):
            for j in range(-dot_size, dot_size + 1):
                dist_to_center = math.sqrt(i ** 2 + j ** 2)
                local_roi.append([xpos + i - 1, ypos + j - 1])
                local_roi_mask.append(dist_to_center < circle_limit)
    return local_roi, local_roi_mask


def check_roi_boundaries(data, values, mask):
    """Pixels of the pencil outside of the map replaced by None."""
    new_values = []
    new_mask = []
    for i in range(len(values)):
        position = values[i]
        if position is not None and position[0] >= 0 and \
                position[1] >= 0 and position[0] < data.nbr_pixels_x and \
                position[1] < data.nbr_pixels_y:
            new_values.append(position)
            new_mask.append(mask[i])
        else:
            new_values.append(None)
            new_mask.append(False)
    return new_values, new_mask


def update_roi_on_multimeshgrids_widget():
    """Update the ROIs on the multimeshgrids widget."""
    # Update the ROIs if displayed on multimeshgrids widget
    if widgets_list.widget_multimeshgrids is not None:
        if shared.exp.display_roi_in_multimeshgrid:
            widgets_list.widget_multimeshgrids.update_MPL()


def ask_user_to_delete_layer_first():
    """Ask the user to delete the glass layer bound to the ROI first."""
    text = "Delete the glass layer bound to the ROI before deleting it."
    msg = QtWidgets.QMessageBox()
    msg.setText("Error : can not delete ROI")
    msg.setInformativeText(text)
    msg.setIcon(QtWidgets.QMessageBox.Critical)
    msg.addButton("OK", QtWidgets.QMessageBox.AcceptRole)
    msg.exec_()
