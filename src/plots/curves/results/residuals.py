from ..plot_curve import PlotCurve
from ....tools import curve_tools

class ResultsResiduals(PlotCurve):
    """Plots the Force - Indentation curve."""

    def __init__(self, parent):
        super().__init__(parent)

        # Labels
        if self.data.curve_distance_units == "nm":
            self.xlabel = "Indentation [nm]"
        elif self.data.curve_distance_units == "um":
            self.xlabel = "Indentation [\u03bcm]"
        self.ylabel = "Force [" + self.data.curve_force_units + "]"

        if self.data.curve_force_units == "nN":
            self.force_factor = 1
        elif self.data.curve_force_units == "pN":
            self.force_factor = 1e3
        if self.data.curve_distance_units == "nm":
            self.x_factor = 1
        elif self.data.curve_distance_units == "um":
            self.x_factor = 1e-3

        if self.data.perform_fit:
            self.get_stiffness_lm_fits()
            self.curve_residuals = self.residuals
            self.segments_fits = None
            self.chisqr = None

        self.define_x_scale()
        self.define_y_scale()
        self.plot_fig()
