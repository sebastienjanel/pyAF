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

"""Widget displaying all the meshgrids."""

import os
import textwrap
from .. import shared
from .. import widgets_list
from PyQt5 import QtGui, QtCore, QtWidgets
from ..tools.gui_tools import PYAFComboBox
from ..tools.gui_tools import PYAFCheckBox
from ..tools.gui_tools import ClearWidgetsFromLayout
from ..tools import misc_tools
from ..tools.PYAFWidget import PYAFWidget
from ..plots.PYAFPlot import PYAFPlot


class MultiMeshgridsWidget(PYAFWidget):
    """Widget which displays all the meshgrids, to be able to compare them.

    If some files are not computed, the meshgrid will be left blanck.
    The meshgrids use the colorscales from the main meshgrids in the results
    widget. If you want to change the color scale for the meshgrids in this
    widget, just change the color scale on the main meshgrid, it will update
    the plots in the multimeshgrid widget.
    """

    def __init__(self, parent):
        super().__init__(
            parent, "widget_multimeshgrids")

        self.current_nbr_cols = 3
        # 10 = small margin
        self.canvas_size = 168.0 + 8

        # Define empty lists to store the elements
        self.qwidget_list = []
        self.canvas_list = []
        self.plots_list = []
        self.title_list = []

        self.VL = QtWidgets.QVBoxLayout()

        self.scroll_area = QtWidgets.QScrollArea()
        self.scroll_area.setWidgetResizable(True)
        self.scroll_area.setMinimumSize(528, 210)  # 528 = 3*176

        # Box with canvases
        self.widget_meshgrids_canvases = QtWidgets.QWidget()
        self.HL_box = QtWidgets.QHBoxLayout()
        self.VL_box = QtWidgets.QVBoxLayout()

        self.grid_layout_box = QtWidgets.QGridLayout()
        self.grid_layout_box.setContentsMargins(0, 0, 0, 0)
        self.grid_layout_box.setHorizontalSpacing(0)
        self.create_meshgrids()

        # Center the gridlayout in the scrollarea
        self.HL_box.addStretch(1)
        self.HL_box.addLayout(self.grid_layout_box)
        self.HL_box.addStretch(1)
        self.VL_box.addStretch(1)
        self.VL_box.addLayout(self.HL_box)
        self.VL_box.addStretch(1)
        self.HL_box.setContentsMargins(0, 0, 0, 0)
        self.VL_box.setContentsMargins(0, 0, 0, 0)
        self.widget_meshgrids_canvases.setLayout(self.VL_box)

        self.scroll_area.setWidget(self.widget_meshgrids_canvases)

        # Box with options
        self.box_options = QtWidgets.QGroupBox("Options")
        HL_box_options = QtWidgets.QHBoxLayout()

        self.list_choser = PYAFComboBox(self, "list_choser")
        self.list_choser.addItem("Piezo height")
        self.list_choser.addItem("Topography")
        self.list_choser.addItem("Elasticity")
        self.list_choser.addItem("Elasticity Slice")
        self.list_choser.addItem("Elasticity (BEC)")
        self.list_choser.addItem("Elasticity Slice (BEC)")
        self.list_choser.addItem("Detachment work")
        self.list_choser.addItem("Detachment force")
        self.list_choser.addItem("Events per curve")
        self.list_choser.addItem("Event max. force")

        self.CB_display_roi = PYAFCheckBox(self, "display_roi", "Display ROIs")

        HL_box_options.addWidget(self.list_choser)
        HL_box_options.addWidget(self.CB_display_roi)
        HL_box_options.addStretch(1)
        self.box_options.setLayout(HL_box_options)

        self.VL.addWidget(self.scroll_area)
        self.VL.addWidget(self.box_options)

        self.setLayout(self.VL)
        self.update_widget()

    def create_meshgrids(self):
        """The plots for the meshgrids are created here and added to the layout."""
        # Clear all widgets from the layout
        ClearWidgetsFromLayout(self.grid_layout_box)

        self.qwidget_list = []
        self.canvas_list = []
        self.plots_list = []
        self.title_list = []

        smallfont = QtWidgets.QApplication.font()
        smallfont.setPointSize(10)

        for datasetid in range(len(shared.exp.list)):
            self.qwidget_list.append(QtWidgets.QWidget())
            # Set the minimum size of the plots, else they can shrink too much
            # when using the resize button of the widget. See bug #358.
            # min_size=[144, 144]
            MPL_canvas = PYAFPlot(
                self,
                "multimeshgrid",
                self.qwidget_list[datasetid],
                [2, 2, 72],
                data_id=datasetid,
                min_size=[144, 144])
            MPL_canvas.setContentsMargins(0, 0, 0, 0)
            MPL_canvas.mpl_connect("scroll_event", self.scroll_on_mesh)
            MPL_canvas.mpl_connect("button_press_event", self.click_on_mesh)

            widg = QtWidgets.QWidget()
            VL = QtWidgets.QVBoxLayout()

            # Get a cleaned up title on two rows
            full_title, title = get_title(datasetid)

            title = QtWidgets.QLabel(title)
            self.title_list.append(title)
            title.setFont(smallfont)
            HL_title = QtWidgets.QHBoxLayout()
            HL_title.addStretch(1)
            HL_title.addWidget(title)
            HL_title.addStretch(1)
            VL.addLayout(HL_title)
            VL.addWidget(MPL_canvas)
            widg.setLayout(VL)

            # For very long titles, add a tooltip
            widg.setToolTip(full_title)

            # Store the canvases for the plot updates
            self.canvas_list.append(MPL_canvas)

            # Store the layouts for the positionning
            self.plots_list.append(widg)

        self.update_positions()

    def update_MPL(self):
        """Update the plots."""
        for datasetid in range(len(shared.exp.list)):
            self.canvas_list[datasetid].update_plot()

    def update_widget(self):
        """Updates the plots and the GUI."""
        # Check first if no data was added or deleted, If it's the case, just
        # recreate the plots and readd them to the layout
        if len(self.qwidget_list) != len(shared.exp.list):
            self.create_meshgrids()

        meshgrid_type = shared.exp.multi_meshgrid_type
        meshgrid_type = misc_tools.get_meshgrid_type_by_id(meshgrid_type)
        self.list_choser.setCurrentIndex(meshgrid_type)

        self.update_MPL()

        self.update_positions()

        # Update the titles
        for datasetid in range(len(shared.exp.list)):
            full_title, title = get_title(datasetid)
            self.title_list[datasetid].setText(title)
            self.plots_list[datasetid].setToolTip(full_title)

        # Display or hide ROIs
        self.CB_display_roi.setChecked(shared.exp.display_roi_in_multimeshgrid)

    def checkbox_clicked(self, name):
        """Called when a checkbox is clicked."""
        if name == "display_roi":
            val = self.CB_display_roi.isChecked()
            shared.exp.display_roi_in_multimeshgrid = val

            self.update_MPL()

    def update_positions(self):
        """Update the positions of the meshgrids depending on the window size."""
        # Remove all widgets from the layout
        for datasetid in range(len(shared.exp.list)):
            self.plots_list[datasetid].setParent(None)

        area_width = self.scroll_area.width()
        nbr_cols = int(area_width / self.canvas_size)

        i = 0
        j = 0
        for datasetid in range(len(shared.exp.list)):
            self.grid_layout_box.addWidget(self.plots_list[datasetid], j, i)
            if i == nbr_cols - 1:
                j = j + 1
                i = 0
            else:
                i = i + 1

        self.current_nbr_cols = nbr_cols

    def scroll_on_mesh(self, event):
        """Fetch scroll events on the meshgrids and scroll all the frame.

        The matplotib plots fetch the scroll events. We have to connect the
        plot to this method to fetch the scroll events. The scroll step is
        multiplied by 120 as it was divided by 120 in matplotlib.
        Then a dummy QWheelEvent is created with only the step value filled;
        this event is then sent to the vertical scrollbar.
        """
        step = event.step * 120.0

        event = QtGui.QWheelEvent(QtCore.QPoint(), step,
                                  QtCore.Qt.NoButton, QtCore.Qt.NoModifier)

        scrollbar = self.scroll_area.verticalScrollBar()
        scrollbar.event(event)

    def click_on_mesh(self, event):
        """Called when clicking on meshgrid.

        Allows to change to another file and to another curve.
        """
        if event.button == 1:
            for plot in self.canvas_list:
                if event.inaxes == plot.canvas.axes:
                    data_id = plot.data_id

                    # Change the id
                    shared.exp.id_selected = data_id

                    xpos, ypos = misc_tools.get_position_on_meshgrid(event)
                    shared.exp.meshgrid_click_xpos = xpos - 1
                    shared.exp.meshgrid_click_ypos = ypos - 1

                    # Access file_changed through widget_compute
                    widgets_list.widget_main.file_changed(
                        option="box", option2="multimeshgrid")

    def list_updated(self, what):
        """Update the list."""
        if what == "list_choser":
            # Change only if another type has been selected
            # (Prevents unwanted refreshing)
            meshgrid_type = misc_tools.get_meshgrid_type_by_id(
                shared.exp.multi_meshgrid_type)

            if self.list_choser.currentIndex() != meshgrid_type:
                val = self.list_choser.currentIndex()
                shared.exp.multi_meshgrid_type = \
                    misc_tools.get_meshgrid_type_as_string(val)

                # Updates
                self.update_MPL()

    def resizeEvent(self, event):
        """Catch the resize event.

        Will refresh the position of the meshgrids if the windows is resized.
        """
        # Unused
        _ = event

        area_width = self.scroll_area.width()
        nbr_cols = int(area_width / self.canvas_size)

        if nbr_cols != self.current_nbr_cols:
            self.update_positions()


def get_title(datasetid):
    """Returns a nicely formated title on two rows, and the full title.

    Uses maximum 30 characters.
    """
    # Wrap title if its too long
    # Allow only 30 characters
    data = shared.exp.list[datasetid]
    title = str(os.path.basename(data.filename))[0:30]
    full_title = title
    title = "\n".join(textwrap.wrap(title, 15))

    # Add an empty line anyway in case of smaller titles
    if len(title) < 15:
        title += "\n"

    return full_title, title
