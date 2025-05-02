# sto_controller.py
import numpy as np
from sto_model import STOModel
from sto_view import STOView

class STOController:
    def __init__(self, view: STOView):
        self.view = view
        self.model = STOModel()
        self.view.btn_calc.clicked.connect(self.update_plot)

    def update_plot(self):
        # 1) Fetch inputs
        weight, thrust = self.view.get_inputs()

        # 2) Create a thrust array for the x-axis
        thr_min = max(1.0, thrust / 10)
        thr_max = thrust * 3
        thrust_array = np.linspace(thr_min, thr_max, 300)

        # 3) Compute STO curves for three weight cases
        sto_center = self.model.compute_STO_curve(weight, thrust_array)
        sto_low    = self.model.compute_STO_curve(weight - 10000, thrust_array)
        sto_high   = self.model.compute_STO_curve(weight + 10000, thrust_array)
        curves = [sto_low, sto_center, sto_high]

        # 4) Compute the highlighted point for the exact inputs
        sto_pt = self.model.compute_STO(weight, thrust)

        # 5) Ask the view to plot everything
        self.view.plot(thrust_array, curves, thrust, sto_pt)