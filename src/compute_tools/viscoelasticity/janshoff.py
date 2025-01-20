"""
Python application for fitting Hertz model to different force curves
"""

import numpy as np
import pandas as pd
from ..stiffness import stiffness_models as sm
from ...tools import fit_tools
from ...tools import curve_tools


class ComputeJanshoff:

    def __init__(self):
        self.coeff = None
        self.vo = None
        self.to = None
        self.tm = None

    def main_compute_janshoff(self, i, j, app, ret, dt, fit_params, o_params):

        tip_geometry = o_params["tip_geometry"]

        if tip_geometry == 0:
            # Sphere
            fit_model_app = self.janshoff_sphere_approach
            fit_model_ret = self.janshoff_sphere_retract
            fit_model_force = sm.hertz_sphere
            fit_model_force_lmfit = sm.hertz_sphere_lmfit
            div = np.sqrt(fit_params["tip_radius"] * 1e-9)
            self.coeff = (4.0 / 3.0) * div / (1.0 - fit_params["poisson_ratio"] ** 2)

        elif tip_geometry == 1:
            # Cone
            fit_model_app = self.janshoff_cone_pyramid_approach
            fit_model_ret = self.janshoff_cone_pyramid_retract
            self.coeff = o_params["coeff"]

        elif tip_geometry == 2:
            # Pyramid
            fit_model_app = self.janshoff_cone_pyramid_approach
            fit_model_ret = self.janshoff_cone_pyramid_retract
            self.coeff = o_params["coeff"]

        elif tip_geometry == 3:
            # 4 Sided Pyramid
            fit_model_app = self.janshoff_cone_pyramid_approach
            fit_model_ret = self.janshoff_cone_pyramid_retract
            div = np.tan(np.radians(fit_params["tip_angle"]))
            self.coeff = 1.342 * ((1.0 - fit_params["poisson_ratio"] ** 2) / div)

        elif tip_geometry == 4:
            # Cylinder
            fit_model_app = self.janshoff_cylinder_approach
            fit_model_ret = self.janshoff_cylinder_retract
            div = 2 * fit_params["tip_radius"] * 1e-9
            self.coeff = (1.0 - fit_params["poisson_ratio"] ** 2) / div

        # Get curve.
        app, ret, _, _, _, smoothing_error = curve_tools.get_curve(None, [i, j], app, ret, mode="compute", dt=dt)

        delta, data_to_fit, final_fit, first, stop = fit_tools.prepare_data_hertzfit(app, dt, fit_params)

        if len(delta) <= 1 or len(data_to_fit) <= 1:
            # Check if there is no data to fit data. Skip if no data is found.
            return

        else:

            # A setpoint that later will be used to select the data for the fit is defined
            setpoint = fit_params["trig_threshold"] * dt["deflection_sensitivity"] * dt["spring_constant"]

            data = {"data_to_fit": data_to_fit[stop::], "delta": delta[stop::]}
            df = pd.DataFrame(data)

            start_f = 0
            stop_f = 33 / 100 * setpoint

            filtered_df = df[(start_f <= df.data_to_fit) & (df.data_to_fit <= stop_f)]

            masked_data_to_fit = np.r_[data_to_fit[0:stop], filtered_df["data_to_fit"].values]
            masked_delta = np.r_[delta[0:stop], filtered_df["delta"].values]

            delta0_fit, e0_fit = fit_tools.do_hertz_fit(fit_model_force, fit_model_force_lmfit,
                                                        masked_delta, masked_data_to_fit)

            poc_index_fit = (np.abs(delta - delta0_fit)).argmin() + first
            real_poc = [app[0][poc_index_fit] * 1e09, app[1][poc_index_fit] * 1e09]

    # Approach segment models ##########################################################################################

    def janshoff_sphere_approach(self, time, e0, beta):
        a = np.sqrt(np.pi) * np.special.gamma(1 - beta)
        b = 4 * np.special.gamma(5/2 - beta)
        y = 3 * e0 * a * (1 / self.coeff * b) * np.power(self.vo * time, 3/2)
        return y

    def janshoff_cone_pyramid_approach(self, time, e0, beta):
        y = e0 * 1 / (self.coeff * (2 - 3 * beta + np.power(beta, 2))) * np.power(self.vo * time, 2)
        return y

    def janshoff_cylinder_approach(self, time, e0, beta):
        y = e0 * 1 / (self.coeff * (1 - beta)) * self.vo * time
        return y

    # Retract segment models ###########################################################################################

    def solved_integral(self, time, beta):
        a = 2 * np.power(time, 2 - beta) / np.power(beta, 2) - 3 * beta + 2
        b = np.power(2, 1 / 1 - beta) * (1 - beta) * 2 * np.power(np.power(time - self.tm, 1 - beta), 1 / 1 - beta) \
            + (beta - 2) * time
        c1 = np.power(2, 1 / 1 - beta)
        c2 = np.power(np.power(time - self.tm, 1 - beta), 1 / 1 - beta)
        c = (1 - beta) * (2 - beta) * np.power(c1 * c2, beta - 1)
        return a - (b / c)

    def janshoff_sphere_retract(self, time, e0, beta):
        y = np.power(self.vo, 3 / 2) * e0 * np.power(self.to, beta) * 1 / self.coeff * self.solved_integral(time, beta)
        return y

    def janshoff_cone_pyramid_retract(self, time, e0, beta):
        y = np.power(self.vo, 2) * e0 * np.power(self.to, beta) * 1 / self.coeff * self.solved_integral(time, beta)
        return y

    def janshoff_cylinder_retract(self, time, e0, beta):
        y = self.vo * e0 * np.power(self.to, beta) * 1 / self.coeff * self.solved_integral(time, beta)
        return y

    # Approach segment models ##########################################################################################

    def janshoff_sphere_approach_lmfit(self, time, e0, beta):
        pass

    def janshoff_cone_pyramid_approach_lmfit(self, time, e0, beta):
        pass

    def janshoff_cylinder_approach_lmfit(self, time, e0, beta):
        pass

    # Retract segment models ###########################################################################################

    def janshoff_sphere_retract_lmfit(self, time, e0, beta):
        pass

    def janshoff_cone_pyramid_retract_lmfit(self, time, e0, beta):
        pass

    def janshoff_cylinder_retract_lmfit(self, time, e0, beta):
        pass