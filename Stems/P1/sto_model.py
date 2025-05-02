# sto_model.py
import numpy as np

class STOModel:
    def __init__(self):
        # Constants for takeoff calculation
        self.S = 1000.0       # Wing area (ft^2)
        self.CL_max = 2.4     # Maximum lift coefficient
        self.CD = 0.0279      # Drag coefficient
        self.rho = 0.002377   # Air density (slugs/ft^3)
        self.gc = 32.174      # Gravitational constant (lbm·ft)/(lbf·s²)

    def compute_V_stall(self, weight):
        return np.sqrt(weight / (0.5 * self.rho * self.S * self.CL_max))

    def compute_V_TO(self, weight):
        return 1.2 * self.compute_V_stall(weight)

    def compute_A(self, thrust, weight):
        return self.gc * thrust / weight

    def compute_B(self, weight):
        return self.gc / weight * (0.5 * self.rho * self.S * self.CD)

    def compute_STO(self, weight, thrust):
        """
        Compute takeoff distance for a single thrust and weight via the analytical integral.
        STO = (1/(2B)) * ln(A / (A - B * V_TO**2))
        """
        Vto = self.compute_V_TO(weight)
        A = self.compute_A(thrust, weight)
        B = self.compute_B(weight)
        if A - B * Vto**2 <= 0:
            return np.nan
        return 1.0/(2*B) * np.log(A / (A - B * Vto**2))

    def compute_STO_curve(self, weight, thrust_array):
        """
        Vectorized STO calculation over an array of thrust values.
        """
        Vto = self.compute_V_TO(weight)
        A = self.gc * thrust_array / weight
        B = self.gc / weight * (0.5 * self.rho * self.S * self.CD)
        inner = A - B * Vto**2
        # Only valid where inner > 0
        sto = np.where(inner > 0, 1.0/(2*B) * np.log(A / inner), np.nan)
        return sto