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

"""Layer for images in VTK."""

import os
from scipy import ndimage
from ... import shared
from .layer import Layer
from ...tools import misc_tools
from ..utils.common import fill_image_sizes
from ..utils.vtk_image import get_actor_plane, image_source, read_image


class LayerSingleImage(Layer):
    """Creates an single image layer.

    Works only with a single tiff file (8 or 16 bits, gray or color).
    Sometimes the tiff you will get will not have the right orientation
    (it will be mirored or inverted ...) Then you will have to specify
    a tiff orientation in the menu
    See http://www.awaresystems.be/imaging/tiff/tifftags/orientation.html
    """

    def __init__(self, parent, path):
        super().__init__()

        self.parent = parent
        self.layer_id = len(shared.layer_list)
        self.type = "single_image"
        self.filename = os.path.basename(path)
        self.filepath = path
        self.bind = True
        self.from_itk = None  # DEBUG value

        # 8 or 16 bits, gray or color (not the stacks)
        images = read_image(path)
        if images is None:
            self.load_error = True
            return None
        self.pixels, self.gray = images

        size = [self.pixels.shape[1], self.pixels.shape[0]]
        # Store the sizes of the images (are defined in layer.py)
        fill_image_sizes(self, size, misc_tools.get_tiff_calibration(path))

        # Textures: the image, or its gray levels through the LUT
        self.image = image_source(self.pixels)
        self.image_gray = image_source(self.gray)

        self.actor = self.get_actor()

        self.init_axes()

        self.addToRenderer()

    def get_actor(self):
        """Plane with the image (gray levels colored by the LUT)."""
        lut = self.get_current_lut_options()
        flip = self.get_current_flip_options()
        if self.smooth_gaussian > 0:
            # Gaussian filter (standard deviation in pixels), each color
            pixels = self.gray if lut["use"] else self.pixels
            sigma = [self.smooth_gaussian] * 2 + [0] * (pixels.ndim - 2)
            image = image_source(ndimage.gaussian_filter(pixels, sigma))
        else:
            image = self.image_gray if lut["use"] else self.image
        return get_actor_plane(image, [self.size_x, self.size_y],
                               self.pos_z, lut, flip)

    def addToRenderer(self):
        """Add actor to renderer."""
        self.parent.renderer.AddActor(self.actor)

    def remove_actors(self):
        """Remove actor from layer."""
        self.parent.renderer.RemoveActor(self.actor)

    def update_actor(self):
        """Update the image layer actor."""
        self.remove_actors()

        self.actor = self.get_actor()

        # Add the actor to the renderer
        self.addToRenderer()

        # Render before updating the axes
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

        # Refresh axes
        self.init_axes()

        self.rotate(None, None, new_angle=self.angle, new_actor=True)
        self.translate(0, 0, 0)

        self.hide_or_display_layer()

        self.set_opacity()

        # Update the rendering
        if shared.VTK_first is False:
            self.parent.interactor.GetRenderWindow().Render()

    def set_opacity(self, value=None):
        """Sets the opacity of the actors."""
        if value is None:
            value = self.opacity

        self.actor.GetProperty().SetOpacity(value / 100.0)

        # Call parent's method set the axes opacity
        super().set_opacity(value)
