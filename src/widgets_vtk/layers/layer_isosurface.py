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

"""Isosurface layer."""

from ... import consts
from ... import shared
from .layer import Layer
from ..utils.vtk_isosurface import get_isosurface

if consts.ALLOW_ITK:
    import itk


class LayerIsoSurface(Layer):
    """Isosurface layer.

    Stack_id is the id of the layer to take the tiff stack from.
    """

    def __init__(self, parent, stack_id, iso_file):
        super().__init__()

        self.parent = parent
        self.layer_id = len(shared.layer_list)
        self.type = "isosurface"
        self.filename = "Isosurface"

        self.original_layer = shared.layer_list[stack_id]
        self.isosurf_orig_layer_id = stack_id  # For loading and saving
        if iso_file is None:
            self.input_image = self.original_layer.input_image
        else:
            # Load from a saved tiff file
            image_type = itk.Image[itk.UC, 3]
            reader = itk.ImageFileReader[image_type].New()
            reader.SetFileName(str(iso_file))
            reader.Update()
            self.input_image = reader.GetOutput()

        self.size_x = self.original_layer.size_x
        self.size_y = self.original_layer.size_y
        self.size_z = self.original_layer.size_z
        self.res_x = self.original_layer.res_x
        self.res_y = self.original_layer.res_y
        self.res_z = self.original_layer.res_z

        color = self.isosurface_color
        img = self.input_image
        flip = [self.original_layer.flip_x, self.original_layer.flip_y, 0]
        self.actor = get_isosurface(img, self.sizes, self.res, color, flip)

        self.init_axes()

        self.addToRenderer()

        self.hide_or_display_layer()

        self.set_opacity()

    def addToRenderer(self):
        """Add actor to the rendered."""
        self.parent.renderer.AddActor(self.actor)

    def remove_actors(self):
        """Remove the actor."""
        self.parent.renderer.RemoveActor(self.actor)

    def update_actor(self):
        """Update the actor."""
        self.remove_actors()

        color = self.isosurface_color
        img = self.input_image
        flip = [self.original_layer.flip_x, self.original_layer.flip_y, 0]
        self.actor = get_isosurface(img, self.sizes, self.res, color, flip)

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
        """Sets the opacity of the actor."""
        if value is None:
            value = self.opacity

        self.actor.GetProperty().SetOpacity(value / 100.0)

        # Call parent's method set the axes opacity
        super().set_opacity(value)
