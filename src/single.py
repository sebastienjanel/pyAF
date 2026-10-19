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

"""
The Data class contains the data corresponding to a single file. Each time a
file is loaded, a Data object is created and saved in the experiment's **list**
variable. The **id_selected** variable in Experiment let's you know which file
is currently displayed in pyAF.

The **dontsave** variable in the Experiment and Data classes contain a list of
variables which will not be saved in .pyaf files. Instead of listing all the
variables to be automatically saved (there are too many), we exclude the few
ones which should not be saved. For example, the arrays with the data are
already saved in the temporary hdf5 file, so they do not need to be saved
again. For further informations take a look at the src/load_and_save/save.py
module.

The **tempfile** variable in the Experiment class contains references to the
current temporary hdf5 file. See the src/tools/temp_file.py for more details.

.. note::
   Variables starting with an underscore (_value) are reimplemented
   in the @property below (with getters and setters). See some intersting
   discussion about getters and setters in python here:
   <http://stackoverflow.com/questions/6618002/\
   python-property-versus-getters-and-setters>`_.
"""

import numpy
from scipy import special
from . import widgets_list
from .tools import results_sorting



# Upper limit of the automatic color scales: percentile of the values (the
# maximum would let a few outliers squeeze the colors of the other pixels,
# the values above are drawn with the saturation color)
SCALE_PERCENTILE = 95
# Automatic color scale of the elasticity: median of the values of the
# sample (substrate excluded) multiplied and divided by 10 ** (n sigma),
# sigma being the robust standard deviation of log10(E) (e.g. 0-240 kPa for
# a cell of 15-40 kPa next to glass, cone model); bottom at 0 if below this
# fraction of the top
ELASTICITY_SIGMAS = 4
ELASTICITY_ZERO_FRACTION = 0.1
# Bottom of a log scale starting at 0: this percentile of the positive
# values
ELASTICITY_PERCENTILES = (2, 98)

# Results of the viscoelasticity fit (index in data.visco_results)
VISCO_RESULTS = ("e0", "fluidity", "tc", "f0", "r2", "height")


def robust_max(values):
    """Upper limit of a color scale: percentile of the finite values (their
    maximum if the percentile is not positive), None without values."""
    values = numpy.asarray(values, dtype=float).ravel()
    values = values[numpy.isfinite(values)]
    if values.size == 0:
        return None
    high = float(numpy.percentile(values, SCALE_PERCENTILE))
    if high <= 0:
        high = float(values.max())
    return high

def robust_range(values, log=False):
    """Limits (low, high) of the automatic color scale of the elasticity,
    None without positive values: median of log10(E) +- ELASTICITY_SIGMAS
    robust standard deviations (1.4826 x median absolute deviation), not
    moved by the tail of the distribution and by a minority of other pixels
    (e.g. a substrate not excluded). On a linear scale, the bottom is 0 if
    it is below ELASTICITY_ZERO_FRACTION of the top.
    """
    values = numpy.asarray(values, dtype=float).ravel()
    values = values[numpy.isfinite(values) & (values > 0)]
    if values.size == 0:
        return None
    logs = numpy.log10(values)
    median = float(numpy.median(logs))
    sigma = 1.4826 * float(numpy.median(numpy.abs(logs - median)))
    if sigma == 0:
        # Uniform values: +- 10 %
        sigma = numpy.log10(1.1) / ELASTICITY_SIGMAS
    low = 10.0 ** (median - ELASTICITY_SIGMAS * sigma)
    high = 10.0 ** (median + ELASTICITY_SIGMAS * sigma)
    if not log and low < ELASTICITY_ZERO_FRACTION * high:
        low = 0.0
    return low, high


class SingleFile:
    """Objet containing the data and options for a single file."""

    def __init__(self, parent, filename, temp_file, unique_id):
        # List of things wich are note saved in the pickle dump
        self.dontsave = ["dontsave",
                         "parent",
                         "temp_file",
                         "datadict",
                         "_rupture_force2",
                         "_events_per_curve",
                         "nbr_points_x",  # For retrocompatibility
                         "nbr_points_y"]  # For retrocompatibility

        # Note, the retrocompatibility lines can be removed for 1.6 release.

        self.parent = parent
        self.temp_file = temp_file
        self.unique_id = unique_id

        # Define empty arrays, will be filled in the update method.
        self._rupture_force2 = None
        self._events_per_curve = None

        # A dictionnary containing references to all the data in the HDF5
        # file. Will be filled below.
        self.datadict = {}

        # ---- SAVED by pickle dump from here on ------------------------------
        self.filename = filename

        # Parameters of the dataset
        self._original_spring_constant = None
        self._spring_constant = None
        self._original_deflection_sensitivity = None
        self._deflection_sensitivity = None
        self._original_temperature = None
        self._temperature = None
        self._used_spring_constant = None

        # Define values to know what has been calculated
        self.stiffness_calculated = False
        self.stiffness_corrected = False
        self.work_and_rupture_force1_calculated = False
        self.loading_rates_calculated = False
        self.events_calculated = False
        self.visco_calculated = False

        # --- Header informations ---------------------------------------------
        self.file_type = None
        # The file type can be one of the following :
        # Nanoscope (Single File)
        # Nanoscope (Force Volume)
        # Nanoscope (Peak Force)
        # JPK (Single File)
        # JPK (Force Map)
        # JPK (QI)

        self.microscope_name = None
        self.version = None
        self.scanner_file = None
        self.real_nbr_pixels = None  # For incomplete JPK files
        self.sens_z_scan = None  # Nanoscope
        self.scan_angle = None
        self.scan_size_x = None  # In nm
        self.scan_size_y = None  # In nm
        self.nbr_pixels_x = None
        self.nbr_pixels_y = None
        self.x_size = None
        self.y_size = None
        # True if the z values come from a height sensor, False if they are
        # calculated from the piezo (None for the old .pyaf files)
        self.zsensor_used = None
        # Note, z_size is not defined as it can be different for each pixel
        # in JPK files. Access this value by using the get_x_dist function
        # in the tools/math_tools module.
        self.nbr_points_per_curve_approach = None  # For example 512 or 1024
        self.nbr_points_per_curve_retraction = None
        self.nbr_points_per_curve_approach_real = None  # Less than 512 or 1024
        self.nbr_points_per_curve_retraction_real = None
        self.nbr_points_per_pause_curve = None
        self.nbr_points_per_modulation_curve = None
        self.image_number_lines = None
        self.number_curves = None
        self.ramp_size = None  # In nm
        self.scan_rate = None  # Speed of scan, in hertz
        self.approach_velocity = None
        self.retraction_velocity = None
        self.trig_threshold = None
        self.frame_direction = None  # Up or down
        self.date = None
        # Can be different than nbr_points_(xy)
        self.image_samps_per_line = None
        self.piezo_image_samps_per_line = None
        self.note = None
        self.z_closed_loop = None
        self.xy_closed_loop = None
        self.retracted_delay = None
        self.extended_delay = None
        self.plane_fit = None  # Nanoscope
        self.realtime_plane_fit = None  # Nanoscope
        self.offline_plane_fit = None  # Nanoscope
        self.x_factor = None  # Nanoscope
        self.matrix_length = None  # Nanoscope
        self.pause_duration = None  # JPK
        # PeakForce frequency (Hz) of the Nanoscope PeakForce files
        self.pft_freq = None
        self.segment_durations = None   # JPK
        self.corrupted_curves = None

        # Color tables id's for the meshgrids (default)
        # First value is the id, second is min, third is max,
        # fourth is middle, fifth is extra if needed
        self.color_opts_piezo = [0, 0, 0, 0, 0]
        self.color_opts_topo = [0, 0, 0, 0, 0]
        self.color_opts_stiffness = [2, 0, 0, 0, 0]
        self.color_opts_work = [12, 0, 0, 0, 0]
        self.color_opts_rupture_force1 = [12, 0, 0, 0, 0]
        self.color_opts_events_per_curve = [6, 0, 0, 0, 0]
        self.color_opts_rupture_force2 = [12, 0, 0, 0, 0]
        # Viscoelasticity: E0 as the elasticity (Viridis), fluidity with
        # Inferno, loss tangent with a diverging scale centered on 1 (equal
        # storage and loss moduli: solid-like below, liquid-like above)
        self.color_opts_visco_e0 = [2, 0, 0, 0, 0]
        self.color_opts_visco_fluidity = [3, 0, 1, 0.5, 0]
        self.color_opts_visco_tan_delta = [6, 0, 2, 1, 0]
        # Storage (Viridis, as E0) and loss (Magma) moduli
        self.color_opts_visco_storage = [2, 0, 0, 0, 0]
        self.color_opts_visco_loss = [4, 0, 0, 0, 0]
        # Quality of the elasticity fit: R^2 from 0.9 (white to red: the
        # bad fits in red, the good ones light), and the pixels with the
        # point of contact fixed (Reds: fixed in red)
        self.color_opts_fit_r2 = [14, 0.9, 1, 0.95, 0]
        self.color_opts_fit_poc_fixed = [10, 0, 1, 0.5, 0]
        self.color_saturation = [1, 1, 1, 1]
        self.color_negative = [0, 0, 0, 1]
        self.color_nan = [0.8, 0.8, 0.8, 1]  # grey
        self.max_piezo = None
        self.max_piezo_display = None
        self.max_topo = None
        self.max_stiffness = None
        self.max_work = None
        self.max_rupture_force1 = None
        self.max_nbr_events = None
        self.max_rupture_force2 = None

        # Nan color in 3D can be transparent (is optional)
        self.nan_color_transparent = False

        # Color names
        names = ["Red", "Blue", "Green", "Yellow", "Orange", "Black", "White"]

        # Meshgrid
        self.discard_corrupted_curves_applied = False
        self.roi_list = []
        self.roi_color_names = names
        # Log color scale for the elasticity maps
        self.stiffness_color_log = False
        self.roi_selected_row = None

        # Profiles
        self.profile_list = []
        self.profile_starts = []
        self.profile_count = 0
        self.profiles_color_names = names

        # Smoothing
        self.sg_smoothing_enabled = False
        self.sg_smoothing_order = 3
        self.sg_smoothing_width = 9
        self.sg_smoothing_uniform = True
        self.used_sg_smoothing_enabled = None
        self.used_sg_smoothing_order = None
        self.used_sg_smoothing_width = None
        self.used_sg_smoothing_uniform = None

        # Fitting parameters for the pocs
        self.fitparam_poc_skip_start = 100  # nm
        self.fitparam_poc_fit_length = 500  # nm
        self.fitparam_poc_refit_option = 500  # nm
        self.fitparam_poc_noise_multiplicator = 5.0
        self.fitparam_poc_refit_times = 2
        # Fitting parameters for the jocs
        self.fitparam_joc_skip_start = 100  # nm
        self.fitparam_joc_fit_length = 500  # nm
        self.fitparam_joc_refit_option = 500  # nm
        self.fitparam_joc_noise_multiplicator = 5.0
        self.fitparam_joc_refit_times = 2
        # Fitting parameters for the jocs (events)
        self.fitparam_events_joc_skip_start = 100  # nm
        self.fitparam_events_joc_fit_length = 500  # nm
        self.fitparam_events_joc_refit_option = 500  # nm
        self.fitparam_events_joc_noise_multiplicator = 5.0
        self.fitparam_events_joc_refit_times = 2

        # Used values
        # Fitting parameters for the pocs
        self.used_fitparam_poc_skip_start = None
        self.used_fitparam_poc_fit_length = None
        self.used_fitparam_poc_refit_option = None
        self.used_fitparam_poc_noise_multiplicator = None
        self.used_fitparam_poc_refit_times = None
        # Fitting parameters for the jocs
        self.used_fitparam_joc_skip_start = None
        self.used_fitparam_joc_fit_length = None
        self.used_fitparam_joc_refit_option = None
        self.used_fitparam_joc_noise_multiplicator = None
        self.used_fitparam_joc_refit_times = None
        # Fitting parameters for the jocs (events)
        self.used_fitparam_events_joc_skip_start = None
        self.used_fitparam_events_joc_fit_length = None
        self.used_fitparam_events_joc_refit_option = None
        self.used_fitparam_events_joc_noise_multiplicator = None
        self.used_fitparam_events_joc_refit_times = None

        # Zooms on plots
        self.zoom_fit_preview_factor = 2
        self.auto_zoom_in = False

        # Stiffness
        self.max_indentation_index = 0

        # Events
        self.events_algorithm = "kernel"  # kernel or msf
        self.events_curve_type = "curve"  # curve or detection
        self.display_fit_events_joc2_preview = True
        self.events_kernel_detection_threshold = 5
        self.events_kernel_adaptive_threshold = 5
        self.events_msf_detection_threshold = 1
        self.events_msf_window_size = 60
        self.events_fit_size = 50
        self.display_fit_event = True
        self.display_fit_event_seg = True

        # Events kernels
        self.kernel_size = 2
        self.adaptive_threshold_option = False
        self.adaptive_smoothing_window = 100
        self.adaptive_smoothing_order = 2

        # Used params
        self.used_events_algorithm = None
        self.used_events_kernel_detection_threshold = None
        self.used_events_kernel_adaptive_threshold = None
        self.used_events_msf_detection_threshold = None
        self.used_events_msf_window_size = None
        self.used_events_fit_size = None
        self.used_kernel_size = None
        self.used_kernel = None
        self.used_adaptive_threshold_option = None
        self.used_adaptive_smoothing_window = None
        self.used_adaptive_smoothing_order = None

        self.events_results_display_annotations = 0
        # 0 = none, 1 = annotations
        self.events_display_results_filter_dist = False
        self.events_results_display_joc = False
        self.events_results_display_fit_joc = False
        # Distance filter (see tools.event_filter_tools)
        self.events_results_filter_dist_min = 0
        self.events_results_filter_dist_max = 0
        self.events_results_filter_dist_keep = False
        self.events_results_filter_slope_min = 0
        self.events_results_filter_slope_max = 0

        # Information for the models
        # 0 Tomography, 1 Fit (perform_fit is derived from it)
        self.stiffness_mode = 1
        self.stiffness_model_selected = 0  # Select Hertz by default
        self.tip_radius = 75
        # Radius of the tip entered for each model (see
        # contact_models.DEFAULT_RADIUS for the others)
        self.tip_radii = {}
        # Half-opening angle, and the angle entered for each model (see
        # contact_models.DEFAULT_ANGLE for the others)
        self.tip_angle = 35
        self.tip_angles = {}
        self.poisson_ratio = 0.5
        self.indentation_start = 0
        self.indentation_stop = 0
        self.indentation_step = 0
        self.force_start = 0
        self.force_stop = 100
        # Fit: no fit if the indentation is smaller than this percentage of
        # the piezo travel after the contact (sample too stiff for the
        # cantilever), and no value if the fit moves the point of contact by
        # more than this percentage of the indentation
        self.fit_min_indentation = 5
        self.fit_poc_window = 25
        self.tomography = 0
        self.lr_coef = 0.33
        #self.mfev = 50
        self.fit_range_type = 0

        # Store the parameters used for the computation
        self.used_deflection_sensitivity = None
        self.used_stiffness_model_selected = None
        # Tomography or Fitting, and fit range type (see
        # used_perform_fit and get_used_fit_range_type)
        self.used_stiffness_mode = None
        self.used_fit_range_type = None
        self.used_tip_radius = None
        self.used_tip_angle = None
        self.used_poisson_ratio = None
        self.used_indentation_start = None
        self.used_force_start = None
        self.used_fit_min_indentation = None
        self.used_fit_poc_window = None
        self.used_indentation_stop = None
        self.used_force_stop = None
        self.used_indentation_step = None
        self.used_tomography = None
        self.used_lr_coef = None

        # Viscoelasticity (Ting fit of the approach and the retraction,
        # power-law rheology, see compute_tools/viscoelasticity)
        # Solver: 0 analytical (constant indentation velocities), 1
        # numerical (measured indentation; always used with the bottom
        # effect correction)
        self.visco_solver = 1
        self.visco_bec = False
        # Part of the curve fitted: indentation > -contact offset [nm]
        self.visco_contact_offset = 100
        self.visco_nbr_points = 300
        # Fluidity: automatic initial value, or the given value
        self.visco_beta_auto = True
        self.visco_beta_init = 0.2
        # Offset of the retraction baseline (viscous drag) corrected
        self.visco_drag_correction = True
        # Frequency of the storage and loss moduli maps [Hz] (computed from
        # E0 and the fluidity, no new fit)
        self.visco_frequency = 1.0
        # Results [i, j, k], k: VISCO_RESULTS
        self.visco_results = None
        self.used_visco_solver = None
        self.used_visco_bec = None
        self.used_visco_contact_offset = None
        self.used_visco_nbr_points = None
        self.used_visco_beta_auto = None
        self.used_visco_beta_init = None
        self.used_visco_drag_correction = None

        # Display options
        self.meshgrid_type = "piezo"
        self.display_curve_type = "defl_ext"

        self.display_trace_retrace = 2
        # 0 = Trace
        # 1 = Retrace
        # 2 = Both

        self.display_poc = False
        self.display_fit_poc = False
        self.display_segments = False
        self.display_fits_stiffness = False
        # Indentation depth for the stiffness meshgrid and the 3D
        self.stiffness_depth_view = 0
        # Depth for the stiffness (for a flat slice at a defined height)
        self.stiffness_slice_depth = 1000  # In nm
        # Display the point of jump of contact
        self.display_joc = False
        self.display_fit_joc = False
        self.display_surface = False
        self.display_force = False
        # Display the joc determined by the events algorithm
        self.display_joc_events = False

        # Force units for the curve plot
        self.curve_force_units = "nN"
        self.curve_distance_units = "nm"
        # Meshgrid units
        self.meshgrid_units = "um"

        # Slices options
        self.slice_position_bottom = 0
        self.slice_position_top = 0
        self.slice_position_left = 0
        self.slice_position_right = 0
        self.slice_view_from_direction = 0
        # 0 = Top, 1 = Bottom, 2 = Left, 3 = Right, 4 = 1st Profile
        # The y min/max values let you define the size of the z scale
        self.slices_z_mode = "auto"
        self.slices_z_min = 0
        self.slices_z_max = 1000
        # Display slice without topography
        self.display_slice_with_topo = True
        # Set aspect ratio of the slices
        self.slice_aspect_ratio = True

        # Slices titles
        self.title_slice_top = "Top view, slice %i"
        self.title_slice_bottom = "Bottom view, slice %i"
        self.title_slice_left = "Left view, slice %i"
        self.title_slice_right = "Right view, slice %i"
        self.title_slice_profile = "Profile 1"

        # Slices options Opengl (will be reset after update is called)
        self.opengl_slice_top = 0
        self.opengl_slice_bottom = 0
        self.opengl_slice_left = 0
        self.opengl_slice_right = 0

        self.opengl_display_colorbar = False
        self.opengl_afm_flat_view = False
        self.opengl_surf_type = "Stiffness"

        # Smoothing of vtk tomographies
        # Type : None, "loop", "butterfly"
        self.vtk_smoothing_type = None
        self.vtk_smoothing_iterations = 1
        self.vtk_decimate_target_reduction = 0.2

        # If the user has flattened the data in pyAF
        self.applied_flatten_order = None
        self.flatten_applied = False

        # Planefit of the height maps with a substrate ROI as zero plane (see
        # tools.planefit): index of the ROI in roi_list
        self.planefit_applied = False
        self.planefit_roi = None

        # Preview values for the tilt limits
        self.tilt_limit_1 = 200
        self.tilt_limit_2 = 600
        self.tilt_applied = False
        self.check_tilt_option = "nbr_points"  # nbr_points or slopes
        # Values applied
        self.used_tilt_limit_1 = None
        self.used_tilt_limit_2 = None
        self.used_tilt_applied = None

        # Preview values for the stretch
        self.stretch_app_lim1 = 100
        self.stretch_app_lim2 = 500
        self.stretch_ret_lim1 = 100
        self.stretch_ret_lim2 = 500
        self.stretch_len_app = 1000
        self.stretch_len_ret = 1000
        self.stretch_applied_app = False
        self.stretch_applied_ret = False
        # Values applied
        self.used_stretch_app_lim1 = 100
        self.used_stretch_app_lim2 = 500
        self.used_stretch_ret_lim1 = 100
        self.used_stretch_ret_lim2 = 500
        self.used_stretch_len_app = 1000
        self.used_stretch_len_ret = 1000
        self.used_stretch_applied_app = False
        self.used_stretch_applied_ret = False

        # List of indentations to be displayed in the histogram tableview
        self.indentation_list = []

        # Indentation multiplicator (currently not used)
        self.height_multiplicator = 1.0

        # Height value inputed by the user for the stiffness correction
        self.user_h = None
        self.noise_floor = None

        # Which roi was used for the stiffness correction
        # 0 = None, 1 = first ROI, etc...
        self.roi_glass_id = 0

        self.discarded_curves = None

    def update_events_per_curve(self):
        """Fills the self._events_per_curve array with values.

        These have to be calculated by results utils each time because they
        depend on the distance limits and on the grades given by the user.
        Therefore self._events_per_curve can not be a static array and needs
        to be refilled each time.
        """
        force_type = "events_per_curve"
        results_util = results_sorting.GetResults(force_type=force_type)
        results_util.load_data(data_id=self.unique_id)

        empty = numpy.zeros((self.nbr_pixels_x, self.nbr_pixels_y), int)
        self._events_per_curve = empty

        name = "events_per_curve"
        for i in range(self.nbr_pixels_x):
            for j in range(self.nbr_pixels_y):
                value = results_util.get_data_from_curve(i, j, name)
                self._events_per_curve[i][j] = value

                widgets_list.widget_progressbar.update()

    def update_rupture_force2(self):
        """Fills the self._rupture_force2 array with values.

        These have to be calculated by results utils each time because they
        depend on the distance limits and on the grades given by the user.
        Therefore self._rupture_force2 can not be a static array and needs
        to be refilled each time.
        """
        # Biggest force event in curve (events)
        name = "events_rupture_force"
        results_util = results_sorting.GetResults(force_type=name)
        results_util.load_data(data_id=self.unique_id)

        empty = numpy.zeros(
            (self.nbr_pixels_x,
             self.nbr_pixels_y),
            numpy.float32)

        # Note : use numpy.float32 here, the other arrays for the meshgrids
        # also contain float32
        self._rupture_force2 = empty

        for i in range(self.nbr_pixels_x):
            for j in range(self.nbr_pixels_y):
                value = results_util.get_data_from_curve(i, j, name)
                self._rupture_force2[i][j] = value

    def get_node(self, name, condition=True, wrap=False):
        """Method used to load the HDF5 nodes to the datadict.

        Either the node is loaded into the datatict if it does not exist,
        or loaded from the datadict. This serves as a small cache for the
        nodes. The reason this implemented like this is that loading a lot
        of nodes inside the class directly slows a lot the loading when
        loading a lot of files (see #474). In this case, loading it to the
        datadict is much more efficient.
        """
        ret = None
        if condition:
            tf = self.temp_file.file
            if wrap:
                # HDF5 VLArrays need to be wrapped with singleToMultiWrapper
                # so they can be accessed as a normal 2D array
                try:
                    ret = self.datadict[name]
                except KeyError:
                    node = tf.get_node("/data/_" + str(self.unique_id) + name)
                    ret = singleToMultiWrapper(node, self.nbr_pixels_y)
                    self.datadict[name] = ret
            else:
                try:
                    ret = self.datadict[name]
                except KeyError:
                    # print("debug")
                    ret = tf.get_node("/data/_" + str(self.unique_id) + name)

        return ret

    @property
    def curves_approach(self):
        """Getter for the approach curve.

        piezo_extension = curves_approach[0] (nm)
        deflection = curves_approach[1] (nm)
        """
        return self.get_node("/curves/approach")

    @property
    def curves_retraction(self):
        """Getter for the retraction curve.

        piezo_extension = curves_retraction[0] (nm)
        deflection = curves_retraction[1] (nm)
        """
        return self.get_node("/curves/retraction")

    @property
    def curves_pause(self):
        """Getter for the pause curve.

        piezo_extension = curves_pause[0] (nm)
        deflection = curves_pause[1] (nm)
        """
        condition = self.nbr_points_per_pause_curve is not None
        return self.get_node("/curves/pause", condition)

    @property
    def curves_modulation(self):
        """Getter for the pause curve.

                piezo_extension = curves_pause[0] (nm)
                deflection = curves_pause[1] (nm)
                """
        condition = self.nbr_points_per_modulation_curve is not None
        return self.get_node("/curves/modulation", condition)

    @property
    def piezo_image(self):
        """Piezo image, last position of the fully extended piezo.

        Expressed in nm.
        """
        return self.get_node("/piezo_image/piezo_image")

    @property
    def approach_positions(self):
        """Array containing the position of corrupted parts of a curve.

        pos = approach_positions[i][j] (value for pixel i, j)
        print pos
        > [32, 512] for example
        This means that the data between the indices 0 and 32 should not
        be taken into account for the computations as they are
        corrupted.
        """
        return self.get_node("/positions/approach_positions")

    @property
    def retraction_positions(self):
        """Array containing the position of corrupted parts of a curve.

        > pos = retraction_positions[i][j] (value for pixel i, j)
        > pos
        > [32, 512] for example
        This means that the data between the indices 0 and 32 should not
        be taken into account for the computations as they are
        corrupted.
        """
        return self.get_node("/positions/retraction_positions")

    @property
    def pocs_indices(self):
        """Indice in the curve of the point of contact.

        See the pocs_real value for the real position of the point of
        contact, here it is just the indice of the point just BEFORE the
        real point of contact which is found by interpolation.

        > print pocs_indices[i][j]
        > 432 for example
        """
        condition = self.stiffness_calculated
        return self.get_node("/results/pocs_indices", condition)

    @property
    def pocs_real(self):
        """Real position of the point of contact found by interpolation.

        For each pixel a list is returned, containing the piezo
        extension and the deflection.

        > print pocs_real[i][j]
        > [453.32, 1.43] (in nm)
        """
        condition = self.stiffness_calculated
        return self.get_node("/results/pocs_real", condition)

    @property
    def fits_poc(self):
        """Parameters used to do the fit for the point of contact.

        Fit is defined by the following equation:
        y = a*x + b

        > print fits_poc[i][j]
        > [a, b, error]

        Error is the value used for the noise. (y+errror or y-error)
        """
        condition = self.stiffness_calculated
        return self.get_node("/results/fits_poc", condition)
    
    @property
    def indentation(self):
        """Parameters used to do the fit for the point of contact.

        Fit is defined by the following equation:
        y = a*x + b

        > print fits_poc[i][j]
        > [a, b, error]

        Error is the value used for the noise. (y+errror or y-error)
        """
        condition = self.stiffness_calculated
        return self.get_node("/results/indentation", condition)

    @property
    def stiffness_array(self):
        """Array containing the stiffness values in Pascal.

        Returns a list of values for each pixel. If the stiffness is
        computed with the Kasas model, the list contains the stiffness
        values for each indentation depth:
        > print stiffness_array[i][j]
        > [23455, 324543, 32425345, 423332]
        """
        condition = self.stiffness_calculated
        return self.get_node("/results/stiffness_array", condition, True)

    @property
    def topography(self):
        """Topography is the z position of the point of contact in nm.

        > topography[i][j]
        > 1252.0
        """
        condition = self.stiffness_calculated
        return self.get_node("/results/topography", condition)

    @property
    def jocs1_indices(self):
        """Indice of the first jump of contact.

        This position is the real point where there is no more contact
        between the tip and the substrate (Joc on the left).

        This is the indice of the point before the real JOC1.

        ::
                                                  -
                                                 -
            ---------------------1              2
                                  ------      -
                                        ------
        """
        condition = self.work_and_rupture_force1_calculated
        return self.get_node("/results/jocs1_indices", condition)

    @property
    def jocs2_indices(self):
        """Indice of the second jump of contact.

        This position is the point where the attractive part of the
        retraction curve starts. It is more on the right.

        This is the indice of the point before the real JOC2.

        ::
                                                  -
                                                 -
            ---------------------1              2
                                  ------      -
                                        ------
        """
        condition = self.work_and_rupture_force1_calculated
        return self.get_node("/results/jocs2_indices", condition)

    @property
    def jocs1_real(self):
        """Real position of the jump of contact 1.

        This position is found by interpolation. See jocs1_indices for
        more details.
        For each pixel a list is returned, containing the piezo
        extension and the deflection.

        > print jocs1_real[i][j]
        [453.32, 1.43] (in nm)
        """
        condition = self.work_and_rupture_force1_calculated
        return self.get_node("/results/jocs1_real", condition)

    @property
    def jocs2_real(self):
        """Real position of the jump of contact 2.

        This position is found by interpolation. See jocs2_indices for
        more details.
        For each pixel a list is returned, containing the piezo
        extension and the deflection.

        > print jocs2_real[i][j]
        [123.76, 0.45] (in nm)
        """
        condition = self.work_and_rupture_force1_calculated
        return self.get_node("/results/jocs2_real", condition)

    @property
    def work(self):
        """Value of the work for each pixel (in J)."""
        condition = self.work_and_rupture_force1_calculated
        return self.get_node("/results/work", condition)

    @property
    def rupture_force1(self):
        """Value of maximal rupture force for each pixel."""
        condition = self.work_and_rupture_force1_calculated
        return self.get_node("/results/rupture_force1", condition)

    @property
    def fits_joc(self):
        """Fit parameters for the jump of contact detection.

        It is for each pixel.
        Fit is defined by the following equation:
        y = a*x + b

        > print fits_joc[i][j]
        [a, b, error]
        """
        condition = self.work_and_rupture_force1_calculated
        return self.get_node("/results/fits_joc", condition)

    @property
    def events_positions_middle(self):
        """Indice of the events for each pixel (middle position)"""
        condition = self.events_calculated
        return self.get_node(
            "/results/events_positions_middle", condition, True)

    @property
    def events_positions_start(self):
        """Indice of the events for each pixel (start position).

        This position can be used to define the beggining of an event
        but is only an estimation of this position.
        """
        condition = self.events_calculated
        return self.get_node(
            "/results/events_positions_start", condition, True)

    @property
    def events_positions_stop(self):
        """Indice of the events for each pixel (stop position).

        This position can be used to define the end of an event but is
        only an estimation of this position.
        """
        condition = self.events_calculated
        return self.get_node(
            "/results/events_positions_stop", condition, True)

    @property
    def events_forces(self):
        """Forces for all events for each pixel."""
        condition = self.events_calculated
        return self.get_node("/results/events_forces", condition, True)

    @property
    def events_slopes(self):
        """Slopes after each event, for each pixel.

        This can be used to sort events by slopes. (For example to sort
        out tethers or protein-protein unbinding events).
        """
        condition = self.events_calculated
        return self.get_node("/results/events_slopes", condition, True)

    @property
    def events_jocs2_indices(self):
        """Indice of the "jump of contact" on the retraction curve.

        Used for the events.
        All the events on the right side of the joc2 are not saved
        during the computation. This position allows also to compute the
        distance between the joc2 and the event.
        It is called events_distance in pyAF.

        ::
                                                  -
                                                 -
            ----------------------              2
                                  ------      -
                                        ------
        """
        condition = self.events_calculated
        return self.get_node("/results/events_jocs2_indices", condition)

    @property
    def events_jocs2_real(self):
        """Real position of the "jump of contact".

        It is the position on the retraction curve.
        This value is the real joc2, and is found by interpolation.
        > print events_jocs2_real[i][j][0]
        > [342.31, 0.456] (in nm)
        """
        condition = self.events_calculated
        return self.get_node("/results/events_jocs2_real", condition)

    @property
    def events_forces_distance(self):
        """Distance between the event and the joc2.

        Expressed in meters, for each pixel.
        """
        condition = self.events_calculated
        return self.get_node(
            "/results/events_forces_distance", condition, True)

    @property
    def events_fits_joc(self):
        """Fit parameters for the jump of contact events detection.

        Fit is defined by the following equation:
        y = a*x + b

        > print events_fits_joc[i][j]
        [a, b, error]
        """
        condition = self.events_calculated
        return self.get_node("/results/events_fits_joc", condition)

    @property
    def stiffness_corrected_array(self):
        """Stiffness values corrected for the bottom effect."""
        condition = self.stiffness_corrected
        return self.get_node("/results/stiffness_corrected", condition, True)

    @property
    def loading_rates(self):
        """Loading rates values for each event and for each pixel."""
        condition = self.loading_rates_calculated
        return self.get_node("/results/loading_rates", condition, True)

    @property
    def perform_fit(self):
        """True if the stiffness is computed by fitting (Fitting selected),
        False for Tomography. Derived from stiffness_mode, so that
        the computation always matches the displayed choice."""
        return self.stiffness_mode == 1

    @perform_fit.setter
    def perform_fit(self, value):
        self.stiffness_mode = 1 if value else 0

    @property
    def used_perform_fit(self):
        """True if the stiffness was computed by fitting.

        Falls back on the current choice for the files saved before this
        value was stored.
        """
        if self.used_stiffness_mode is None:
            return self.perform_fit
        return self.used_stiffness_mode == 1

    def get_used_fit_range_type(self):
        """Returns the fit range type used for the stiffness computation.

        Falls back on the current value for the files saved before this
        value was stored.
        """
        if self.used_fit_range_type is None:
            return self.fit_range_type
        return self.used_fit_range_type

    def get_used_indentation(self):
        """Returns the (start, step, tomography) used for the stiffness
        computation, to describe the computed segments.

        Falls back on the current parameters for the files saved before these
        values were stored.
        """
        start = self.used_indentation_start
        step = self.used_indentation_step
        tomography = self.used_tomography
        return (self.indentation_start if start is None else start,
                self.indentation_step if step is None else step,
                self.tomography if tomography is None else tomography)

    def reset_discarded(self):
        """Method used to reset the discarded curve array."""
        empty = numpy.zeros((self.nbr_pixels_x, self.nbr_pixels_y))
        self.discarded_curves = empty

    def update(self, color=True):
        """Updates the color options and resets the datadict.

        The color arguments let's you decide if you want to reset the
        max value of the colorscales. This update is needed after a new
        computation.
        """
        # Reset nodes list, so that it will be repopulated
        self.datadict = {}

        if self.nbr_pixels_x is not None:
            self.opengl_slice_top = self.nbr_pixels_x - 1
            self.opengl_slice_right = self.nbr_pixels_y - 1

        # Piezo image
        # Piezo height: autoscale (0 to the maximum of the file)
        self.max_piezo = numpy.amax(self.piezo_image)
        if color:
            self.color_opts_piezo[2] = self.max_piezo
            self.color_opts_piezo[3] = self.max_piezo / 2

        if self.stiffness_calculated:
            # Update max_stiffness
            array_stiffness = []
            for i in range(self.nbr_pixels_x):
                for j in range(self.nbr_pixels_y):
                    for k in range(len(self.stiffness_array[i][j])):
                        array_stiffness.append(self.stiffness_array[i][j][k])

            self.max_stiffness = numpy.amax(array_stiffness)
            self.max_topo = numpy.nanmax(self.topography)
            if color:
                self.color_opts_topo[2] = self.max_topo
                self.color_opts_topo[3] = self.max_topo / 2.0
                limits = robust_range(self.sample_elasticity_values(),
                                      self.stiffness_color_log)
                if limits is not None:
                    self.color_opts_stiffness[1] = limits[0]
                    self.color_opts_stiffness[2] = limits[1]
                    self.color_opts_stiffness[3] = sum(limits) / 2.0

        # Jump of contact and work
        if self.work_and_rupture_force1_calculated:
            # Update max_work and rupture_force1
            self.max_work = numpy.amax(self.work)
            self.max_rupture_force1 = numpy.amax(self.rupture_force1)
            if color:
                self.color_opts_work[2] = self.max_work
                self.color_opts_work[3] = self.max_work / 2.0
                self.color_opts_rupture_force1[2] = self.max_rupture_force1
                self.color_opts_rupture_force1[3] = \
                    self.max_rupture_force1 / 2.0

        # Viscoelasticity
        if self.visco_calculated and color:
            # Moduli: same automatic scale as the elasticity
            for name in ("visco_e0", "visco_storage", "visco_loss"):
                limits = robust_range(getattr(self, name),
                                      self.stiffness_color_log)
                if limits is not None:
                    opts = getattr(self, "color_opts_" + name)
                    opts[1], opts[2] = limits
                    opts[3] = sum(limits) / 2.0

        # Events
        if self.events_calculated:
            # Note that here, unlike the conditions above, the
            # events arrays have already been updated by the load script

            # Update max_nbr_events
            array = []
            array_events_forces = []
            for i in range(self.nbr_pixels_x):
                for j in range(self.nbr_pixels_y):
                    array.append(len(self.events_forces[i][j]))
                    array_events_forces.append(self._rupture_force2[i][j])

            self.max_nbr_events = numpy.amax(array)
            self.max_rupture_force2 = numpy.amax(array_events_forces)
            if color:
                self.color_opts_events_per_curve[2] = self.max_nbr_events
                self.color_opts_events_per_curve[3] = \
                    int(self.max_nbr_events / 2.0)
                self.color_opts_rupture_force2[2] = self.max_rupture_force2
                self.color_opts_rupture_force2[3] = \
                    self.max_rupture_force2 / 2.0

    @property
    def nbr_points_x(self):
        """Old nbr_points_x getter.

        To be removed in 1.5 version
        """
        print(
            "DeprecationWarning: nbr_points_x is deprecated. "
            "Please use data.nbr_pixels_x.")

        return self.nbr_pixels_x

    @property
    def nbr_points_y(self):
        """Old nbr_points_y getter.

        To be removed in 1.5 version
        """
        print(
            "DeprecationWarning: nbr_points_y is deprecated. "
            "Please use data.nbr_pixels_y.")

        return self.nbr_pixels_y

    @property
    def array_curves_approach(self):
        """Old array_curves_approach getter.

        To be removed in 1.5 version
        """
        print(
            "DeprecationWarning: array_curves_approach is deprecated. "
            "Please use data.curves_approach.")

        return self.curves_approach

    @property
    def array_curves_retraction(self):
        """Old array_curves_retraction getter.

        To be removed in 1.5 version
        """
        print(
            "DeprecationWarning: array_curves_retraction is deprecated. "
            "Please use data.curves_retraction.")

        return self.curves_retraction

    @property
    def array_curves_pause(self):
        """Old array_curves_pause getter.

        To be removed in 1.5 version
        """
        print(
            "DeprecationWarning: array_curves_pause is deprecated. "
            "Please use data.curves_pause.")

        return self.curves_pause

    @property
    def calculated(self):
        """Checks if at least one computation has been done on the dataset."""
        if self.stiffness_calculated or self.events_calculated or \
                self.work_and_rupture_force1_calculated:
            return True

        else:
            return False

    @property
    def is_single(self):
        """Getter to check if the data is a force map or a single curve.

        If it is a single curve, will return True
        """
        ft = self.file_type

        if ft == "Nanoscope (Single File)" or ft == "JPK (Single File)":
            return True
        else:
            return False

    @property
    def spring_constant(self):
        """Returns the current spring constant.

        If the user has manually changed the spring constant of the cantilever
        we use this value (original value is never modified).
        1e-9 is for N/m to N/nm conversion

        The spring constant is returned in N/nm but stored in N/m
        """
        if self._spring_constant is None:
            return self._original_spring_constant * 1e-9
        else:
            return self._spring_constant * 1e-9

    @spring_constant.setter
    def spring_constant(self, value):
        """Setter for the spring constant."""
        self._spring_constant = value

    @property
    def used_spring_constant(self):
        """Returns the used spring constant.

        The spring constant is returned in N/nm but stored in N/m
        """
        if self._used_spring_constant is not None:
            return self._used_spring_constant * 1e-9
        else:
            return None

    @used_spring_constant.setter
    def used_spring_constant(self, value):
        """Sets the used spring constant."""
        self._used_spring_constant = value

    @property
    def deflection_sensitivity(self):
        """Getter for the deflection sensitivity."""
        if self._deflection_sensitivity is None:
            return self._original_deflection_sensitivity
        else:
            return self._deflection_sensitivity

    @deflection_sensitivity.setter
    def deflection_sensitivity(self, value):
        """Setter for the deflection sensitivity."""
        self._deflection_sensitivity = value

    @property
    def temperature(self):
        """Getter for the temperature."""
        default_temp = 21.0 + 273.15

        if self._original_temperature == -1 and self._temperature is None:
            temp = default_temp
        elif self._original_temperature != -1 and self._temperature is None:
            temp = self._original_temperature
        elif self._temperature is not None:
            temp = self._temperature

        return temp

    @temperature.setter
    def temperature(self, value):
        """Setter for the temperature."""
        self._temperature = value

    @property
    def events_per_curve(self):
        """Getter for the events_per_curve array."""
        return self._events_per_curve

    @property
    def rupture_force2(self):
        """Getter for the rupture_force2 array."""
        return self._rupture_force2

    @property
    def pocs_corrected(self):
        """Getter for the topography array. Deprecated.

        To be removed for 1.6 release.
        """
        print (
            "DeprecationWarning: pocs_corrected is deprecated. "
            "Please use data.topography.")
        if self.stiffness_calculated:
            return self.topography
        else:
            return None

    @property
    def force_samples_per_curve(self):
        """Getter for the number of points per curve. Deprecated.

        To be removed for 1.6 release.
        """
        print (
            "DeprecationWarning: force_samples_per_curve is deprecated. "
            "Please use data.nbr_points_per_curve_approach or "
            "data.nbr_points_per_curve_retraction.")
        return self.nbr_points_per_curve_approach

    @property
    def force_samples_per_curve_real(self):
        """Getter for the real number of points per curve. Deprecated.

        To be removed for 1.6 release.
        """
        print (
            "DeprecationWarning: force_samples_per_curve_real is deprecated. "
            "Please use data.nbr_points_per_curve_approach_real or "
            "data.nbr_points_per_curve_retraction_real.")
        return self.nbr_points_per_curve_approach_real

    @property
    def force_samples_per_pause(self):
        """Getter for the force_samples_per_pause. Deprecated.

        To be removed for 1.6 release.
        """
        print (
            "DeprecationWarning: force_samples_per_pause is deprecated. "
            "Please use data.nbr_points_per_pause_curve.")
        return self.nbr_points_per_pause_curve

    def get_visco_result(self, name):
        """Map of a result of the viscoelasticity fit (nan if not
        computed)."""
        if not self.visco_calculated or self.visco_results is None:
            return numpy.full((self.nbr_pixels_x, self.nbr_pixels_y),
                              numpy.nan)
        return self.visco_results[:, :, VISCO_RESULTS.index(name)]

    @property
    def visco_e0(self):
        """Modulus E0 of the power law (at t0 = 1 s) [Pa]."""
        return self.get_visco_result("e0")

    @property
    def visco_fluidity(self):
        """Fluidity (exponent of the power law, 0: solid, 1: liquid)."""
        return self.get_visco_result("fluidity")

    @property
    def visco_tan_delta(self):
        """Loss tangent of the power law: tan(pi * fluidity / 2)."""
        return numpy.tan(numpy.pi * self.visco_fluidity / 2)

    def visco_modulus(self):
        """Complex modulus of the power law at visco_frequency [Pa]:
        E*(w) = E0 Gamma(1 - b) (i w t0) ** b (t0 = 1 s), derived from E0
        and the fluidity b (no measure at this frequency)."""
        beta = self.visco_fluidity
        omega = 2 * numpy.pi * getattr(self, "visco_frequency", 1.0)
        return self.visco_e0 * special.gamma(1 - beta) * omega ** beta

    @property
    def visco_storage(self):
        """Storage modulus E' at visco_frequency [Pa]."""
        return self.visco_modulus() * \
            numpy.cos(numpy.pi * self.visco_fluidity / 2)

    @property
    def visco_loss(self):
        """Loss modulus E'' at visco_frequency [Pa]."""
        return self.visco_modulus() * \
            numpy.sin(numpy.pi * self.visco_fluidity / 2)

    def sample_elasticity_values(self, corrected=False):
        """Elasticity values of the sample, for the automatic color scale:
        without the substrate, i.e. the pixels where the point of contact
        of the fit was fixed (sample too stiff for the cantilever, failed
        fit) and the pixels of the substrate ROIs. All the values if no
        pixel is left."""
        array = self.stiffness_corrected_array if corrected \
            else self.stiffness_array
        exclude = numpy.zeros((self.nbr_pixels_x, self.nbr_pixels_y), bool)
        if self.has_fit_quality():
            exclude |= self.fit_poc_fixed == 1
        for roi in self.roi_list:
            if roi.glass_coeffs is not None:
                for i, j in roi.values:
                    exclude[i, j] = True
        values = [v for i in range(self.nbr_pixels_x)
                  for j in range(self.nbr_pixels_y) if not exclude[i, j]
                  for v in array[i][j]]
        if not numpy.any(numpy.isfinite(values)):
            values = [v for i in range(self.nbr_pixels_x)
                      for j in range(self.nbr_pixels_y) for v in array[i][j]]
        return values

    def get_fit_quality(self, index):
        """Map of the quality of the elasticity fit: R^2 (index 0), or 1 if
        the point of contact was fixed and 0 if fitted (index 1); nan if
        not computed (Tomography, files computed before)."""
        values = numpy.full((self.nbr_pixels_x, self.nbr_pixels_y),
                            numpy.nan)
        if self.stiffness_calculated and self.has_fit_quality():
            values = numpy.asarray(
                self.get_node("/results/fit_quality")[:, :, index], float)
        return values

    def has_fit_quality(self):
        """True if the quality of the elasticity fit is stored (Fitting,
        computed with this version)."""
        if not self.stiffness_calculated or not self.used_perform_fit:
            return False
        return self.temp_file.file.__contains__(
            "/data/_" + str(self.unique_id) + "/results/fit_quality")

    @property
    def fit_r2(self):
        """R^2 of the elasticity fit (fitted points after the contact)."""
        return self.get_fit_quality(0)

    @property
    def fit_poc_fixed(self):
        """1 where the point of contact of the elasticity fit was fixed
        (fallback), 0 where it was fitted."""
        return self.get_fit_quality(1)

    def get_color_opts(self, meshgrid_type=None):
        """Color options (list) of a meshgrid type (the current one by
        default)."""
        mt = self.meshgrid_type if meshgrid_type is None else meshgrid_type
        if mt in ("stiffness", "stiffness_slice", "stiffness_corr",
                  "stiffness_corr_slice"):
            return self.color_opts_stiffness
        names = {"piezo": "color_opts_piezo",
                 "topo": "color_opts_topo",
                 "work": "color_opts_work",
                 "rupture_force": "color_opts_rupture_force1",
                 "events_per_curve": "color_opts_events_per_curve",
                 "events_rupture_force": "color_opts_rupture_force2",
                 "visco_e0": "color_opts_visco_e0",
                 "visco_fluidity": "color_opts_visco_fluidity",
                 "visco_tan_delta": "color_opts_visco_tan_delta",
                 "visco_storage": "color_opts_visco_storage",
                 "visco_loss": "color_opts_visco_loss",
                 "fit_r2": "color_opts_fit_r2",
                 "fit_poc_fixed": "color_opts_fit_poc_fixed"}
        return getattr(self, names[mt])

    @property
    def colortableid(self):
        """Gets the id of the current colortable."""
        return self.get_color_opts()[0]

    @colortableid.setter
    def colortableid(self, value):
        """Sets the id of the current colortable."""
        self.get_color_opts()[0] = value

    @property
    def colortable_min_value(self):
        """Gets the min value of the current colortable."""
        return self.get_color_opts()[1]

    @colortable_min_value.setter
    def colortable_min_value(self, value):
        """Sets the min value of the current colortable."""
        self.get_color_opts()[1] = value

    @property
    def colortable_max_value(self):
        """Gets the max value of the current colortable."""
        return self.get_color_opts()[2]

    @colortable_max_value.setter
    def colortable_max_value(self, value):
        """Sets the max value of the current colortable."""
        self.get_color_opts()[2] = value

    @property
    def colortable_middle_value(self):
        """Gets the middle value of the current colortable."""
        return self.get_color_opts()[3]

    @colortable_middle_value.setter
    def colortable_middle_value(self, value):
        """Sets the middle value of the current colortable."""
        self.get_color_opts()[3] = value

    @property
    def fit_params_poc(self):
        """Convenient getter for the fit parameters."""
        fit_params = {
            "poc_skip_start": self.fitparam_poc_skip_start,
            "poc_fit_length": self.fitparam_poc_fit_length,
            "poc_refit_option": self.fitparam_poc_refit_option,
            "poc_noise_multiplicator": self.fitparam_poc_noise_multiplicator,
            "poc_refit_times": self.fitparam_poc_refit_times}

        return fit_params

    @property
    def fit_params_joc(self):
        """Convenient getter for the fit parameters."""
        fit_params = {
            "joc_skip_start": self.fitparam_joc_skip_start,
            "joc_fit_length": self.fitparam_joc_fit_length,
            "joc_refit_option": self.fitparam_joc_refit_option,
            "joc_noise_multiplicator": self.fitparam_joc_noise_multiplicator,
            "joc_refit_times": self.fitparam_joc_refit_times}

        return fit_params

    @property
    def fit_params_joc_events(self):
        """Convenient getter for the fit parameters."""
        fit_params = {
            "joc_skip_start": self.fitparam_events_joc_skip_start,
            "joc_fit_length": self.fitparam_events_joc_fit_length,
            "joc_refit_option": self.fitparam_events_joc_refit_option,
            "joc_noise_multiplicator":
            self.fitparam_events_joc_noise_multiplicator,
            "joc_refit_times": self.fitparam_events_joc_refit_times}

        return fit_params

    @property
    def used_fit_params_poc(self):
        """Convenient getter for the fit parameters."""
        fit_params = {
            "poc_skip_start": self.used_fitparam_poc_skip_start,
            "poc_fit_length": self.used_fitparam_poc_fit_length,
            "poc_refit_option": self.used_fitparam_poc_refit_option,
            "poc_noise_multiplicator":
            self.used_fitparam_poc_noise_multiplicator,
            "poc_refit_times": self.used_fitparam_poc_refit_times}

        return fit_params

    @property
    def used_fit_params_joc(self):
        """Convenient getter for the fit parameters."""
        fit_params = {
            "joc_skip_start": self.used_fitparam_joc_skip_start,
            "joc_fit_length": self.used_fitparam_joc_fit_length,
            "joc_refit_option": self.used_fitparam_joc_refit_option,
            "joc_noise_multiplicator":
            self.used_fitparam_joc_noise_multiplicator,
            "joc_refit_times": self.used_fitparam_joc_refit_times}

        return fit_params

    @property
    def used_fit_params_joc_events(self):
        """Convenient getter for the fit parameters."""
        fit_params = {
            "joc_skip_start": self.used_fitparam_events_joc_skip_start,
            "joc_fit_length": self.used_fitparam_events_joc_fit_length,
            "joc_refit_option": self.used_fitparam_events_joc_refit_option,
            "joc_noise_multiplicator":
            self.used_fitparam_events_joc_noise_multiplicator,
            "joc_refit_times": self.used_fitparam_events_joc_refit_times}

        return fit_params


class singleToMultiWrapper:
    """Class used for the wrapping of one dimensional arrays to two dimensions."""

    def __init__(self, data, nbr_pixels_y):
        """Init as an emtpy array and fill with data."""
        self.array = []
        self.nbr_pixels_y = nbr_pixels_y

        for row in data.iterrows():
            self.array.append(row)

    def __len__(self):
        """Define the __len__, the returned array should have this."""
        return self.nbr_pixels_y

    def __getitem__(self, key):
        """Return the row, depending on the position."""
        return (
            self.array[key * self.nbr_pixels_y:(key + 1) * self.nbr_pixels_y]
        )
