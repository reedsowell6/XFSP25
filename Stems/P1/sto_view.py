# sto_view.py
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

class STOView(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("Takeoff Distance Calculator")
        main_layout = QtWidgets.QVBoxLayout()

        form = QtWidgets.QFormLayout()
        self.le_weight = QtWidgets.QLineEdit()
        self.le_weight.setText("56000")
        self.le_thrust = QtWidgets.QLineEdit()
        self.le_thrust.setText("13000")
        form.addRow("Weight (lb):", self.le_weight)
        form.addRow("Thrust (lb):", self.le_thrust)

        self.btn_calc = QtWidgets.QPushButton("Calculate")
        form.addRow(self.btn_calc)

        main_layout.addLayout(form)

        # Matplotlib canvas
        self.figure = plt.Figure()
        self.canvas = FigureCanvas(self.figure)
        main_layout.addWidget(self.canvas)

        self.setLayout(main_layout)

    def get_inputs(self):
        weight = float(self.le_weight.text())
        thrust = float(self.le_thrust.text())
        return weight, thrust

    def plot(self, thrust_array, curves, thrust_point, sto_point):
        """
        Draw the three STO vs thrust curves and mark the selected point.
        """
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        labels = ["Weight - 10k", "Weight", "Weight + 10k"]
        for curve, label in zip(curves, labels):
            ax.plot(thrust_array, curve, label=label)

        # mark the selected operating point
        ax.plot(thrust_point, sto_point, 'o', markersize=8, label="Selected Point")

        ax.set_xlabel("Thrust (lbf)")
        ax.set_ylabel("Takeoff Distance (ft)")
        ax.legend()
        ax.grid(True)
        self.canvas.draw()