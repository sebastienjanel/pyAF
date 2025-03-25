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

"""Widget used to display slices in the tomography."""

import os
from PyQt5 import QtCore, QtWidgets
from ..plots.PYAFPlot import PYAFPlot
from ..widgets_small.slices_axes import SlicesAxesWidget
from ..tools.gui_tools import PYAFCheckBox
from ..tools.PYAFWidget import PYAFWidget
from ..widgets.progressbar import Progressbar
from .. import widgets_list
from .. import shared
import tables


class SlicesWidget(PYAFWidget):
    """This widget displays the stiffness tomographies."""

    def __init__(self, parent):
        super().__init__(parent, "widget_slices")

        self.canvas_resolution = self.parent.canvas_resolution
        self.canvas_size = [432, 432]
        sizes = [self.canvas_size[0] / self.canvas_resolution,
                 self.canvas_size[1] / self.canvas_resolution,
                 self.canvas_resolution]

        self.data = shared.exp.current_data

        self.VL = QtWidgets.QVBoxLayout()

        # Box with canvases
        self.box_slices = QtWidgets.QGroupBox("Slices")
        self.VL_box_slices = QtWidgets.QVBoxLayout()

        self.W_canvas = QtWidgets.QWidget()
        self.W_canvas.setFixedSize(self.canvas_size[0], self.canvas_size[1])
        self.canvas = PYAFPlot(self, "meshgrid_slices", self.W_canvas, sizes)
        self.canvas.mpl_connect("button_press_event", self.canvas_press)

        # Checkbox to display slice without the topo (flat surface)
        self.CB_display_topo = PYAFCheckBox(self, "topo", "Use topography")

        # Checkbox to set the aspect of the plot
        self.CB_aspect = PYAFCheckBox(self, "aspect", "Force aspect ratio")

        # Lists to chose the slice options
        self.HL_lists = QtWidgets.QHBoxLayout()
        self.label_view_from = QtWidgets.QLabel()
        self.label_view_from.setText("View from")
        self.label_slice_position = QtWidgets.QLabel()
        self.label_slice_position.setText("Slice")
        self.list_slice = QtWidgets.QComboBox()
        self.list_slice.addItem("Top")
        self.list_slice.addItem("Bottom")
        self.list_slice.addItem("Left")
        self.list_slice.addItem("Right")
        self.list_slice.addItem("1st Profile")
        self.list_slice.addItem("1st Profile (Inverted)")
        self.list_slice.activated.connect(
            lambda: self.list_updated("list_slice_view_direction"))
        self.list_slice.setCurrentIndex(self.data.slice_view_from_direction)
        self.list_position = QtWidgets.QComboBox()
        self.list_position.setSizeAdjustPolicy(
            QtWidgets.QComboBox.AdjustToContents)
        self.list_position.activated.connect(
            lambda: self.list_updated("list_position"))
        self.HL_lists.addWidget(self.label_view_from)
        self.HL_lists.addWidget(self.list_slice)
        self.HL_lists.addWidget(self.label_slice_position)
        self.HL_lists.addWidget(self.list_position)
        self.HL_lists.addStretch(1)

        self.GL_slices = QtWidgets.QGridLayout()
        self.GL_slices.addWidget(self.canvas, 0, 0, 1, 1)
        self.GL_slices.addWidget(self.CB_display_topo, 1, 0)
        self.GL_slices.addWidget(self.CB_aspect, 1, 1)
        self.GL_slices.addLayout(self.HL_lists, 2, 0, 1, 0)

        self.VL_box_slices.addLayout(self.GL_slices)
        self.box_slices.setLayout(self.VL_box_slices)

        self.VL.addWidget(self.box_slices)
        self.VL.addStretch(1)
        self.setLayout(self.VL)

        self.update_MPL()

        self.update_widget()

    def closeEvent(self, event):
        """Called when the widget is closed.

        The profile plotted on the meshgrid in the results tab needs to be
        removed.
    """

        # Update the profile on the parent (will remove it)
        try:
            # Try to remove it, if the app is quited from this window
            # MPL_canvas1 will no more exist ...
            self.parent.MPL_canvas1.canvas.update_blit("all", closing=True)

        except tables.exceptions.ClosedNodeError:
            # The refresh will happen but the hdf5 file is close, resulting in
            # a tables.exceptions.ClosedNodeError error. Do nothing in this
            # case.

            pass

        super().closeEvent(event)

    def update_widget(self):
        """Updates the GUI."""
        # Slice (Top, Bottom, Left, Right)
        self.list_slice.setCurrentIndex(self.data.slice_view_from_direction)

        # Update the number of slices in the drop down list
        # Clear the list
        self.list_position.clear()

        direct = self.data.slice_view_from_direction
        if direct == 0 or direct == 1:
            length = self.data.nbr_pixels_y + 1
        elif direct == 2 or direct == 3:
            length = self.data.nbr_pixels_x + 1
        else:
            length = None

        if length is not None:
            for i in range(1, length):
                self.list_position.addItem(str(i))

        ls = self.list_position
        if direct == 0:
            ls.setCurrentIndex(self.data.slice_position_top)
        elif direct == 1:
            ls.setCurrentIndex(self.data.slice_position_bottom)
        elif direct == 2:
            ls.setCurrentIndex(self.data.slice_position_left)
        elif direct == 3:
            ls.setCurrentIndex(self.data.slice_position_right)

        self.CB_display_topo.setChecked(self.data.display_slice_with_topo)
        self.CB_aspect.setChecked(self.data.slice_aspect_ratio)

    def display_options_widget(self):
        """Display the widget with the color options."""
        self.parent.button_clicked("button_meshgrid_display_options")

    def open_single_figure(self):
        """Open the figure in a separate window."""
        PYAFPlot(self, "meshgrid_slices")

    def set_title(self):
        """Method to set a title for the plot."""
        value, ok = QtWidgets.QInputDialog.getText(self, "Set title", "New value",
                                               text=self.data.title_slice_top)

        if ok:
            direct = self.data.slice_view_from_direction

            if direct == 0:
                self.data.title_slice_top = value
            elif direct == 1:
                self.data.title_slice_bottom = value
            elif direct == 2:
                self.data.title_slice_left = value
            elif direct == 3:
                self.data.title_slice_right = value
            elif direct == 4 or direct == 5:
                self.data.title_slice_profile = value

            # Update the plot
            self.update_MPL()

    def save_slices_as(self, mode=None):
        """Save all the slices in a folder."""
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(self,
            self.tr("Choose a file name"), ".", self.tr("TIFF (*.tif)"))

        if not filename:
            return
        else:
            # Add extension (for Linux, on OS X it will work out of the box)
            ends1 = filename.endswith(".tif")
            ends2 = filename.endswith(".tiff")
            if ends1 and ends2 is False:
                filename += ".tif"

            orig_save_path = str(filename).split(".")[0]
            ext = str(filename).split(".")[1]

        # Save the oldposition to reset it afterwards
        old_pos = self.list_position.currentIndex()

        # Create a progressbar
        Progressbar()
        widgets_list.widget_progressbar.set_label("Saving ...")
        widgets_list.widget_progressbar.set_range(0, self.data.nbr_pixels_x)

        # Loop through the slices
        direct = self.data.slice_view_from_direction
        if direct == 0 or direct == 1:
            nbr_slices = self.data.nbr_pixels_x
        elif direct == 2 or direct == 3:
            nbr_slices = self.data.nbr_pixels_y

        for i in range(nbr_slices):
            if direct == 0:
                self.data.slice_position_top = i
            elif direct == 1:
                self.data.slice_position_bottom = i
            elif direct == 2:
                self.data.slice_position_left = i
            elif direct == 3:
                self.data.slice_position_right = i

            save_path = orig_save_path + "_" + str(i) + "." + ext
            label = "Saving " + os.path.basename(save_path) + " ..."
            widgets_list.widget_progressbar.set_label(label)

            PYAFPlot(self, "meshgrid_stiffness_slices", save_path=save_path)

            if mode == "meshgrid":
                save_path = orig_save_path + "_meshgrid_" + str(i) + "." + ext
                PYAFPlot(self, "meshgrid", save_path=save_path)

            widgets_list.widget_progressbar.update()

        # Reset the position
        if self.data.slice_view_from_direction == 0:
            self.data.slice_position_top = old_pos
        elif self.data.slice_view_from_direction == 1:
            self.data.slice_position_bottom = old_pos
        elif self.data.slice_view_from_direction == 2:
            self.data.slice_position_left = old_pos
        elif self.data.slice_view_from_direction == 3:
            self.data.slice_position_right = old_pos

        # Close the progressbar
        widgets_list.widget_progressbar.close()

    def update_MPL(self):
        """Update the plots."""
        self.canvas.update_plot()
        # Update the profile on the parent's meshgrid
        self.parent.MPL_canvas1.canvas.update_blit("all")

    def list_updated(self, what):
        """Called when a list is updated."""
        if what == "list_position":
            direct = self.data.slice_view_from_direction
            index = self.list_position.currentIndex()

            if direct == 0:
                self.data.slice_position_top = index
            elif direct == 1:
                self.data.slice_position_bottom = index
            elif direct == 2:
                self.data.slice_position_left = index
            elif direct == 3:
                self.data.slice_position_right = index

            # Update the plot
            self.update_MPL()

        elif what == "list_slice_view_direction":
            direct = self.list_slice.currentIndex()
            self.data.slice_view_from_direction = direct

            if direct == 4 or direct == 5:
                # Check for profile
                if self.data.profile_list == []:
                    QtWidgets.QMessageBox.warning(self, "Error",
                                              "Create a profile first !")
                    return False
                # Profile slice, no need to change the slice number
                self.list_position.setEnabled(False)
            else:
                self.list_position.setEnabled(True)

            # Update the plots
            self.update_MPL()

        self.update_widget()

    def canvas_press(self, event):
        """Called when there is a click on the canvas.

        Will open a menu on right click, or mouve the position on the
        currently selected pixel on left click.
        """
        if event.button == 1 and event.inaxes == self.canvas.canvas.axes:
            # Left click, select the pixel (curve) on the meshgrid
            xpos, ypos = self.get_position_on_meshgrid(event)
            if xpos is not None:
                widgets_list.widget_main.change_curve(xpos, ypos)

        elif event.button == 3 and event.inaxes == self.canvas.canvas.axes:
            # Right click, open a menu
            # Get position in screen coordinates (top left position of canvas)
            globalpos = self.canvas.mapToGlobal(QtCore.QPoint(0, 0))
            # Recalculate positions with canvas position
            x = globalpos.x() + event.x
            y = globalpos.y() + (self.canvas_size[1] - event.y)
            pos = QtCore.QPoint(x, y)
            self.popUpMenu(pos)

    def get_position_on_meshgrid(self, event):
        """Get the position on the meshgrid when clicking on the canvas.

        The position is returned in user coordinates (1 to X + 1).
        """
        xsize = self.data.x_size
        ysize = self.data.y_size

        if self.data.slice_view_from_direction == 0:
            return int((event.xdata * 1000.0) / xsize) + 1, \
                self.data.slice_position_top + 1
        elif self.data.slice_view_from_direction == 1:
            return int((event.xdata * 1000.0) / xsize) + 1, \
                self.data.slice_position_bottom + 1
        elif self.data.slice_view_from_direction == 2:
            return self.data.slice_position_left + 1, \
                int((event.xdata * 1000.0) / ysize) + 1
        elif self.data.slice_view_from_direction == 3:
            return self.data.slice_position_right + 1, \
                int((event.xdata * 1000.0) / ysize) + 1
        else:
            # Disabled for personal profiles, but could be implemented
            return None, None

    def set_z_scale_size(self):
        """Display a widget allowing to change the z scale."""
        axes_widget = SlicesAxesWidget(self)
        axes_widget.setWindowTitle("Set z scale")
        axes_widget.resize(200, 100)
        axes_widget.activateWindow()
        axes_widget.show()

    def hide_or_display_red_line(self):
        """Hide or display the red profile lines on the plot."""
        if shared.exp.display_line_stiffness_slice_at_height:
            shared.exp.display_line_stiffness_slice_at_height = False
        else:
            shared.exp.display_line_stiffness_slice_at_height = True
        self.update_MPL()

    def popUpMenu(self, pos):
        """Method to open a pop up menu.

        Called with right click on figure.
        """
        # Create menu
        menu = QtWidgets.QMenu()

        # Define actions
        action = QtWidgets.QAction("Open figure", self)
        action.triggered.connect(self.open_single_figure)
        menu.addAction(action)

        action = QtWidgets.QAction("Color options", self)
        action.triggered.connect(self.display_options_widget)
        menu.addAction(action)

        if shared.exp.display_line_stiffness_slice_at_height:
            text = "Hide red line"
        else:
            text = "Show red line"
        action = QtWidgets.QAction(text, self)
        action.triggered.connect(self.hide_or_display_red_line)
        menu.addAction(action)

        action = QtWidgets.QAction("Set z scale size", self)
        action.triggered.connect(self.set_z_scale_size)
        menu.addAction(action)

        action = QtWidgets.QAction("Set title", self)
        action.triggered.connect(self.set_title)
        menu.addAction(action)

        # Add action to save all figures
        if (self.data.slice_view_from_direction != 4
                or self.data.slice_view_from_direction != 5):
            menu.addSeparator()
            action = QtWidgets.QAction("Save slices as ...", self)
            action.triggered.connect(self.save_slices_as)
            menu.addAction(action)
            action = QtWidgets.QAction("Save slices (+ meshgrid) as ...", self)
            action.triggered.connect(lambda: self.save_slices_as("meshgrid"))
            menu.addAction(action)

        # Display the menu
        menu.popup(pos, menu.menuAction())
        menu.exec_()

    def checkbox_clicked(self, name):
        """Called when a checkbox is clicked."""
        if name == "topo":
            val = self.CB_display_topo.isChecked()
            self.data.display_slice_with_topo = val
            self.update_MPL()
        elif name == "aspect":
            self.data.slice_aspect_ratio = self.CB_aspect.isChecked()
            self.update_MPL()
