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

"""Check of an installation or of an application built with PyInstaller
(see pyAF.spec): run pyAF with the environment variable PYAF_SMOKE_TEST=1.
The resources, the theme, the libraries of the computations and of the 3D
rendering are loaded, then pyAF quits with 0 if everything works.
PYAF_SMOKE_FILE: a file to load, whose elasticity is computed (processes of
the computations); PYAF_SMOKE_STACK: an image stack added to the 3D scene;
PYAF_SMOKE_RENDER: image of the 3D scene (needs OpenGL); PYAF_SMOKE_3D and
PYAF_SMOKE_3D_DIR: images loaded in the real 3D window (needs a screen);
PYAF_SMOKE_LOG: file where the results are also written."""

import os
import tempfile
import traceback


def run(app):
    """Runs the checks, prints them and returns the exit code (0: all the
    checks passed)."""
    failures = []
    # Also written in the file PYAF_SMOKE_LOG (application without console
    # on Windows: no output)
    log = os.environ.get("PYAF_SMOKE_LOG")

    def output(text):
        print(text, flush=True)
        if log:
            with open(log, "a", encoding="utf-8") as file:
                file.write(text + "\n")

    def check(name, function):
        try:
            detail = function()
            output("OK   " + name + (" (" + detail + ")" if detail else ""))
        except Exception:
            failures.append(name)
            output("FAIL " + name + "\n" + traceback.format_exc())

    def resources():
        from . import misc_tools
        from ..gui_styles.theme_handler import get_icons_path
        paths = [misc_tools.get_resource_path(name) for name in
                 ("fix.png", "icon.png", "01.png", "triangle_down.png",
                  "dark.qss", "light.qss", "dark_plot.mplstyle")]
        paths += [misc_tools.get_font_path(),
                  os.path.join(get_icons_path(), "check.svg")]
        missing = [path for path in paths if not os.path.isfile(path)]
        if missing:
            raise FileNotFoundError(", ".join(missing))
        return str(len(paths)) + " files"

    def theme():
        from ..gui_styles.theme_handler import theme_handler
        theme_handler("Dark")
        theme_handler("Light")
        if not app.styleSheet():
            raise RuntimeError("no stylesheet")

    def plot():
        import matplotlib
        matplotlib.use("QtAgg")
        from matplotlib import font_manager
        from matplotlib.figure import Figure
        from . import misc_tools
        font_manager.fontManager.addfont(misc_tools.get_font_path())
        figure = Figure()
        figure.add_subplot(111).plot([0, 1], [0, 1])
        figure.canvas.draw()
        return "matplotlib " + matplotlib.__version__

    def science():
        import numpy
        import scipy.optimize
        import lmfit
        import pandas
        import statsmodels.api  # noqa: F401
        import seaborn  # noqa: F401
        import tifffile  # noqa: F401
        result = scipy.optimize.curve_fit(lambda x, a: a * x,
                                          numpy.arange(5.0),
                                          2 * numpy.arange(5.0))[0]
        if abs(result[0] - 2) > 1e-6:
            raise RuntimeError("wrong fit")
        return "numpy " + numpy.__version__ + ", lmfit " + \
            lmfit.__version__ + ", pandas " + pandas.__version__

    def hdf5():
        import numpy
        import tables
        path = os.path.join(tempfile.mkdtemp(prefix="pyaf_smoke_"), "t.h5")
        with tables.open_file(path, "w") as file:
            file.create_carray("/", "a", tables.Float32Atom(), (10, 10),
                               filters=tables.Filters(complevel=9,
                                                      complib="blosc"))
            file.root.a[:] = numpy.ones((10, 10))
        with tables.open_file(path) as file:
            if float(file.root.a[:].sum()) != 100:
                raise RuntimeError("wrong values")
        return "tables " + tables.__version__

    def vtk():
        # The classes of the 3D view, and their implementations: without
        # them, nothing is drawn, without error at the import
        from vtkmodules.vtkCommonCore import vtkVersion
        from vtkmodules.vtkRenderingCore import vtkTextRenderer
        from ..widgets_vtk import vtk_api
        for name, expected in ((vtk_api.vtkPolyDataMapper(), "OpenGL"),
                               (vtk_api.vtkRenderer(), "OpenGL"),
                               (vtkTextRenderer(), "FreeType")):
            if expected not in name.GetClassName():
                raise RuntimeError(name.GetClassName() + ": not " + expected)
        return "VTK " + vtkVersion.GetVTKVersion()

    def scene():
        # Actors of the 3D view, from an image stack
        import numpy
        from ..widgets_vtk import vtk_api
        from ..widgets_vtk.utils import common, vtk_image
        from ..widgets_vtk.utils import vtk_isosurface
        renderer = vtk_api.vtkRenderer()
        z, y, x = numpy.mgrid[-16:16, -16:16, -16:16]
        stack = 255 * (x ** 2 + y ** 2 + z ** 2 < 12 ** 2).astype(numpy.uint8)
        iso = vtk_isosurface.get_isosurface(stack, [32.0, 32.0], [1.0] * 3,
                                            [0.9, 0.3, 0.2], [1, 0, 0])
        renderer.AddActor(iso)
        points = iso.GetMapper().GetInput().GetNumberOfPoints()
        if not points:
            raise RuntimeError("empty isosurface")
        detail = str(points) + " points, "
        pixels = numpy.tile(numpy.arange(64, dtype=numpy.uint8) * 4, (64, 1))
        plane = vtk_image.get_actor_plane(
            vtk_image.image_source(pixels), [48.0, 48.0], -14.0,
            {"use": True, "min_color": [0.0, 0.0, 0.0],
             "max_color": [0.0, 1.0, 0.0], "threshold": 0, "invert": False},
            {"flip_x": True, "flip_y": False})
        renderer.AddActor(plane)
        # A real image stack, as in the stack layer, one plane per frame
        stack_path = os.environ.get("PYAF_SMOKE_STACK")
        if stack_path:
            images = vtk_image.read_frames(stack_path)
            if images is None:
                raise RuntimeError("stack not read: " + stack_path)
            for index, frame in enumerate(images[0]):
                renderer.AddActor(vtk_image.get_actor_plane(
                    vtk_image.image_source(frame), [40.0, 40.0], index * 2.0,
                    {"use": False}, {"flip_x": False, "flip_y": False}))
            detail += "stack " + str(images[0].shape) + ", "
        lut = vtk_api.vtkLookupTable()
        lut.Build()
        bar = vtk_api.vtkScalarBarActor()
        bar.SetLookupTable(lut)
        common.use_app_font(bar.GetLabelTextProperty())
        renderer.AddActor(bar)
        axes = vtk_api.vtkCubeAxesActor()
        axes.SetBounds(plane.GetBounds())
        axes.SetCamera(renderer.GetActiveCamera())
        renderer.AddActor(axes)
        detail += "3D scene"
        # Rendering: needs OpenGL (not on the servers of the builds), the
        # image can be compared between two builds
        path = os.environ.get("PYAF_SMOKE_RENDER")
        if path:
            from vtkmodules.vtkRenderingCore import vtkRenderWindow
            window = vtkRenderWindow()
            window.SetOffScreenRendering(1)
            window.SetSize(600, 400)
            window.SetMultiSamples(0)
            window.AddRenderer(renderer)
            renderer.ResetCamera()
            window.Render()
            shot = vtk_api.vtkWindowToImageFilter()
            shot.SetInput(window)
            shot.Update()
            writer = vtk_api.vtkPNGWriter()
            writer.SetFileName(path)
            writer.SetInputConnection(shot.GetOutputPort())
            writer.Write()
            detail += ", rendered in " + path
        return detail

    def filters():
        # Filters of the stacks: values computed with ITK, whose filters
        # they replace (voxels of the isosurface at the level 1)
        import numpy
        from ..widgets_vtk.utils import stack_filters
        z, y, x = numpy.mgrid[-16:16, -16:16, -16:16]
        stack = (numpy.mod(x * x + 2 * y * y + 3 * z * z, 251) *
                 (x * x + y * y + z * z < 14 ** 2)).astype(numpy.uint8)
        binary = stack_filters.binary_threshold(stack, 50, 255, 0, 255)
        if int(binary.sum()) != 2333250:
            raise RuntimeError("wrong threshold")
        for sigma, voxels, total in ((1.0, 18402, 2307123),
                                     (2.0, 24520, 2289064)):
            smooth = stack_filters.discrete_gaussian(binary, sigma)
            if int((smooth >= 1).sum()) != voxels or \
                    abs(int(smooth.sum()) / total - 1) > 1e-3:
                raise RuntimeError("wrong smoothing for sigma " + str(sigma))

    def tiff():
        # Reading of the images of the 3D view (8 and 16 bits, gray and
        # color, one or several frames)
        import numpy
        import tifffile
        from ..widgets_vtk.utils import vtk_image
        folder = tempfile.mkdtemp(prefix="pyaf_smoke_")
        ramp = (numpy.arange(3 * 16 * 20).reshape(3, 16, 20) % 251).astype(
            numpy.uint8)
        cases = {"gray8": (ramp, "minisblack", ramp.shape),
                 "gray16": (ramp.astype(numpy.uint16) * 100, "minisblack",
                            ramp.shape),
                 "rgb": (numpy.stack([ramp // 2, ramp, ramp // 3], axis=-1),
                         "rgb", ramp.shape + (3,))}
        for name, (values, photometric, shape) in cases.items():
            path = os.path.join(folder, name + ".tif")
            tifffile.imwrite(path, values, photometric=photometric)
            images = vtk_image.read_frames(path)
            if images is None:
                raise RuntimeError(name + ": not read")
            if images[0].shape != shape or images[1].shape != ramp.shape:
                raise RuntimeError(name + ": shape " + str(images[0].shape))
            # Gray levels: the values, from the min to the max for 16 bits
            expected = ramp if name != "gray16" else \
                (ramp.astype(float) * 255.0 / ramp.max()).astype(numpy.uint8)
            if not numpy.array_equal(images[1], expected):
                raise RuntimeError(name + ": wrong values")
        return str(len(cases)) + " formats"

    def interface():
        from .. import consts, pyaf  # noqa: F401
        if not consts.ALLOW_VTK:
            raise RuntimeError("3D view not available")
        from ..compute_tools import microrheology  # noqa: F401
        from ..compute_tools.stiffness import stiffness_fit  # noqa: F401

    def wait(milliseconds):
        """Runs the events of the interface during this time."""
        from qtpy import QtCore
        loop = QtCore.QEventLoop()
        QtCore.QTimer.singleShot(milliseconds, loop.quit)
        loop.exec()

    def view_3d(main_window):
        # Real 3D window (needs a screen): the images of PYAF_SMOKE_3D
        # (type=path, separated by the path separator), and images of the
        # window in PYAF_SMOKE_3D_DIR
        images = os.environ.get("PYAF_SMOKE_3D")
        if not images:
            return
        import numpy
        from PIL import Image
        from ..widgets_vtk import vtk_api
        from ..widgets_vtk.main import OpenglWidget
        from ..widgets_vtk.utils.create_layer import create_layer
        folder = os.environ["PYAF_SMOKE_3D_DIR"]
        window = OpenglWidget(main_window, [0])
        wait(1500)
        for number, item in enumerate(images.split(os.pathsep)):
            kind, path = item.split("=", 1)
            create_layer(kind, path)
            canvas = window.canvas
            # As the user sees it (no reset of the camera), unless asked
            if os.environ.get("PYAF_SMOKE_3D_RESET"):
                canvas.renderer.ResetCamera()
            canvas.interactor.GetRenderWindow().Render()
            wait(1500)
            shot = vtk_api.vtkWindowToImageFilter()
            shot.SetInput(canvas.interactor.GetRenderWindow())
            shot.Update()
            writer = vtk_api.vtkPNGWriter()
            name = os.path.join(folder, str(number) + "_" + kind + ".png")
            writer.SetFileName(name)
            writer.SetInputConnection(shot.GetOutputPort())
            writer.Write()
            pixels = numpy.asarray(Image.open(name).convert("L"), float)
            output(kind + " " + os.path.basename(path) + ": " +
                   str(canvas.renderer.GetActors().GetNumberOfItems()) +
                   " actors, image " + str(pixels.shape) + " mean " +
                   format(pixels.mean(), ".1f") + " max " +
                   format(pixels.max(), ".0f"))
        window.close()

    def projects():
        # Saved projects: classes saved under the module names of src and of
        # pyAF.src, the one of the running code in the applications
        from . import pickle_tools
        from .. import experiment
        for module in ("src.experiment", "pyAF.src.experiment"):
            data = ("c" + module + "\nResultGroup\n.").encode("ascii")
            if pickle_tools.loads(data) is not experiment.ResultGroup:
                raise RuntimeError("wrong class for " + module)
        return "module names"

    def computation():
        # File given with PYAF_SMOKE_FILE: loaded, and its elasticity
        # computed with the processes of the computations (they re-run the
        # application built with PyInstaller)
        import numpy
        from qtpy import QtCore, QtWidgets
        from .. import pyaf, shared, widgets_list
        from ..load_and_save.pre_check import CheckFiles
        path = os.environ["PYAF_SMOKE_FILE"]
        # No dialog (no user): the message boxes are accepted
        QtWidgets.QMessageBox.exec = lambda box: QtWidgets.QMessageBox.Ok
        QtWidgets.QMessageBox.exec_ = QtWidgets.QMessageBox.exec
        settings = QtCore.QSettings()
        old_temp = settings.value("tempPath")
        settings.setValue("tempPath",
                          tempfile.mkdtemp(prefix="pyaf_smoke_") + os.sep)
        try:
            window = pyaf.MainWindow(CheckFiles([path]).infos)
            window.show()
            wait(500)
            widgets_list.widget_compute.BT_compute_stiffness.click()
            wait(500)
            data = shared.exp.list[0]
            values = [data.stiffness_array[i][j][0]
                      for i in range(data.nbr_pixels_x)
                      for j in range(data.nbr_pixels_y)
                      if len(data.stiffness_array[i][j])]
            view_3d(window)
            window.close()
        finally:
            if old_temp is None:
                settings.remove("tempPath")
            else:
                settings.setValue("tempPath", old_temp)
        if not values:
            raise RuntimeError("no elasticity value")
        return str(len(values)) + " curves, " + str(shared.exp.cores) + \
            " processes, median E " + \
            format(float(numpy.median(values)) / 1e3, ".3g") + " kPa"

    checks = [("resources", resources), ("theme", theme), ("plot", plot),
              ("science", science), ("hdf5", hdf5), ("vtk", vtk),
              ("filters", filters), ("tiff", tiff), ("scene", scene),
              ("projects", projects), ("interface", interface)]
    if os.environ.get("PYAF_SMOKE_FILE"):
        checks.append(("computation", computation))
    for name, function in checks:
        check(name, function)

    output("SMOKE TEST " + ("FAILED: " + ", ".join(failures) if failures
                            else "PASSED"))
    return 1 if failures else 0
