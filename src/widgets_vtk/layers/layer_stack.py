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

"""Tiff stack layer."""

import os
import itk
import vtk
from ... import shared
from .layer import Layer
from ..utils.common import fill_image_sizes
from ...tools import misc_tools
from ..utils.vtk_image import get_actor_plane, image_source, read_frames


class LayerStack(Layer):
    """Tiff stack layer."""

    def __init__(self, parent, path):
        super().__init__()

        self.parent = parent
        self.layer_id = len(shared.layer_list)
        self.type = "stack"
        self.filename = os.path.basename(path)
        self.filepath = path
        self.bind = True

        # 8 or 16 bits, gray or color: pixels to display, and their gray
        # levels (for the LUT, the filters and the isosurfaces, as an ITK
        # image)
        images = read_frames(path)
        if images is None:
            # The image was not recognized
            self.load_error = True
            return None
        self.pixels, self.gray = images
        self.input_image_orig = itk.image_from_array(self.gray)
        self.input_image = self.input_image_orig

        frames, rows, cols = self.gray.shape
        # Store the sizes of the images (are defined in layer.py)
        fill_image_sizes(self, [cols, rows, frames],
                         misc_tools.get_tiff_calibration(path))

        self.build_actors()

        self.init_axes()

        self.addToRenderer()

    def build_actors(self):
        """One plane per image of the stack: the pixels, or their gray
        levels if filtered or through the LUT."""
        # Apply filters to the image if needed
        self.filter_image()

        lut = self.get_current_lut_options()
        filtered = self.use_binary_thresholding or self.smooth_gaussian > 0
        if lut["use"] or filtered:
            images = itk.array_from_image(self.input_image)
        else:
            images = self.pixels

        flip = self.get_current_flip_options()
        size = [self.size_x, self.size_y]
        self.actor = vtk.vtkAssembly()
        self.actor_list = []
        for i, image in enumerate(images):
            actor = get_actor_plane(image_source(image), size,
                                    i * self.res_z, lut, flip)
            self.actor_list.append(actor)
            self.actor.AddPart(actor)

    def addToRenderer(self):
        """Add actor to renderer."""
        self.parent.renderer.AddActor(self.actor)

    def remove_actors(self):
        """Remove the actors."""
        self.parent.renderer.RemoveActor(self.actor)

    def update_actor(self):
        """Update the actor."""
        self.remove_actors()
        self.build_actors()

        self.init_axes()
        self.addToRenderer()

        self.rotate(None, None, new_angle=self.angle, new_actor=True)
        self.translate(0, 0, 0)

        self.hide_or_display_layer()

        self.set_opacity()

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def set_opacity(self, value=None):
        """Sets the opacity of the actors.

        You can set a global opacity for all the images in the stack, or
        define a specific opacity for one single image (and set another
        opacity for the other images)
        """
        for i in range(len(self.actor_list)):
            if i != self.slide_pos:
                val = self.others_opacity / 100.0
                self.actor_list[i].GetProperty().SetOpacity(val)

        val = self.single_opacity / 100.0
        self.actor_list[self.slide_pos].GetProperty().SetOpacity(val)

        # Call parent's method set the axes opacity
        super().set_opacity(value)

    def filter_image(self):
        """Applies the filters (binary threshold, then smoothing) to the
        image."""
        it = itk.Image[itk.UC, 3]
        self.input_image = self.input_image_orig

        if self.use_binary_thresholding:
            filter_type = itk.BinaryThresholdImageFilter[it, it]
            threshold_filter = filter_type.New()

            threshold_filter.SetInput(self.input_image)
            threshold_filter.SetLowerThreshold(self.lower_thresh)
            threshold_filter.SetUpperThreshold(self.upper_thresh)
            threshold_filter.SetOutsideValue(self.lower_thresh_value)
            threshold_filter.SetInsideValue(self.upper_thresh_value)

            self.input_image = threshold_filter.GetOutput()
            self.input_image.Update()

        if self.smooth_gaussian > 0:
            filter_type = itk.DiscreteGaussianImageFilter[it, it]
            smooth_filter = filter_type.New()

            smooth_filter.SetInput(self.input_image)
            # Standard deviation in pixels
            smooth_filter.SetVariance(self.smooth_gaussian ** 2)

            self.input_image = smooth_filter.GetOutput()
            self.input_image.Update()
