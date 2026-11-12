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

"""The part of VTK used by pyAF, imported as vtk by the 3D view.

"import vtk" loads all the libraries of VTK, which the application built
with PyInstaller then contains (about 450 MB, most of them not used). A class
used by the 3D view has to be added here.
"""

# noqa: F401 for the whole list: names for the other modules
from vtkmodules.vtkCommonCore import (  # noqa: F401
    VTK_MAJOR_VERSION, VTK_UNSIGNED_CHAR, VTK_FONT_FILE, mutable, vtkPoints,
    vtkUnsignedCharArray, vtkLookupTable, vtkStringArray, vtkDoubleArray)
from vtkmodules.vtkCommonDataModel import (  # noqa: F401
    vtkTriangle, vtkPolyData, vtkCellArray, vtkImageData, vtkPolygon,
    vtkCellLocator)
from vtkmodules.vtkCommonExecutionModel import vtkTrivialProducer  # noqa: F401
from vtkmodules.vtkCommonTransforms import vtkTransform  # noqa: F401
from vtkmodules.vtkFiltersCore import (  # noqa: F401
    vtkCleanPolyData, vtkTubeFilter, vtkPolyDataNormals,
    vtkPolyDataConnectivityFilter, vtkMarchingCubes, vtkDecimatePro,
    vtkAppendPolyData)
from vtkmodules.vtkFiltersGeneral import vtkTransformFilter  # noqa: F401
from vtkmodules.vtkFiltersModeling import (  # noqa: F401
    vtkLoopSubdivisionFilter, vtkFillHolesFilter)
from vtkmodules.vtkFiltersSources import vtkSphereSource, vtkPlaneSource  # noqa: F401
from vtkmodules.vtkFiltersTexture import vtkTextureMapToPlane  # noqa: F401
from vtkmodules.vtkIOImage import vtkPNGWriter  # noqa: F401
from vtkmodules.vtkImagingCore import vtkImageFlip, vtkImageMapToColors  # noqa: F401
from vtkmodules.vtkInteractionStyle import (  # noqa: F401
    vtkInteractorStyleTrackballCamera)
from vtkmodules.vtkRenderingAnnotation import (  # noqa: F401
    vtkScalarBarActor, vtkCubeAxesActor)
from vtkmodules.vtkRenderingCore import (  # noqa: F401
    vtkPolyDataMapper, vtkActor, vtkTexture, vtkWindowToImageFilter,
    vtkRenderer, vtkPropPicker, vtkProperty, vtkLightKit, vtkAssembly)
# Implementations of the abstract classes above, registered when imported
# (else, nothing is drawn, without any error at the import): OpenGL for the
# rendering, FreeType for the text
import vtkmodules.vtkRenderingFreeType  # noqa: F401
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401
