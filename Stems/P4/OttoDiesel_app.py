# OttoDiesel_app.py
#region imports
from OttoDiesel_GUI import Ui_Form
from PyQt5 import QtWidgets as qtw
import sys
from Otto import ottoCycleController
from Diesel import dieselCycleController
from Dual import dualCycleController    # <-- new
from Air import *
# Matplotlib integration
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
#endregion

class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.calculated = False

        # create canvas for plotting
        self.figure = Figure(figsize=(8,8), tight_layout=True, frameon=True, facecolor='none')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot()
        self.main_VerticalLayout.addWidget(self.canvas)

        # connect signals
        self.rdo_Metric.toggled.connect(self.setUnits)
        self.btn_Calculate.clicked.connect(self.calcCycle)
        self.cmb_Abcissa.currentIndexChanged.connect(self.doPlot)
        self.cmb_Ordinate.currentIndexChanged.connect(self.doPlot)
        self.chk_LogAbcissa.stateChanged.connect(self.doPlot)
        self.chk_LogOrdinate.stateChanged.connect(self.doPlot)
        self.cmb_OttoDiesel.currentIndexChanged.connect(self.selectCycle)

        # instantiate controllers
        self.otto   = ottoCycleController()
        self.diesel = dieselCycleController()
        self.dual   = dualCycleController()
        self.controller = self.otto

        # collect widgets
        self.someWidgets = []
        self.someWidgets += [self.lbl_THigh, self.lbl_TLow, self.lbl_P0, self.lbl_V0, self.lbl_CR]
        self.someWidgets += [self.le_THigh, self.le_TLow, self.le_P0, self.le_V0, self.le_CR]
        # dual-cycle inputs (only if UI generated them)
        if hasattr(self, 'lbl_PR') and hasattr(self, 'le_PR') and hasattr(self, 'lbl_CO') and hasattr(self, 'le_CO'):
            self.someWidgets += [self.lbl_PR, self.le_PR, self.lbl_CO, self.le_CO]
        # state display
        self.someWidgets += [self.le_T1, self.le_T2, self.le_T3, self.le_T4,
                              self.lbl_T1Units, self.lbl_T2Units, self.lbl_T3Units, self.lbl_T4Units]
        # performance readouts
        self.someWidgets += [self.le_PowerStroke, self.le_CompressionStroke,
                              self.le_HeatAdded, self.le_Efficiency,
                              self.lbl_PowerStrokeUnits, self.lbl_CompressionStrokeUnits,
                              self.lbl_HeatInUnits]
        # plotting controls
        self.someWidgets += [self.rdo_Metric, self.cmb_Abcissa, self.cmb_Ordinate,
                              self.chk_LogAbcissa, self.chk_LogOrdinate,
                              self.ax, self.canvas]

        # hand off widgets to controllers
        self.otto.setWidgets(self.someWidgets)
        self.diesel.setWidgets(self.someWidgets)
        # only set dual widgets if UI defines them (matches guard above)
        if hasattr(self, 'lbl_PR') and hasattr(self, 'le_PR') and hasattr(self, 'lbl_CO') and hasattr(self, 'le_CO'):
            self.dual.setWidgets(self.someWidgets)

        self.show()
        # initial populate
        self.controller.updateView()

    def doPlot(self):
        self.controller.updateView()

    def selectCycle(self):
        choice = self.cmb_OttoDiesel.currentText().lower()
        if 'otto' in choice:
            self.controller = self.otto
        elif 'diesel' in choice:
            self.controller = self.diesel
        else:
            self.controller = self.dual
        self.gb_Input.setTitle(f'Input for Air Standard {self.controller.model.name}')
        self.controller.updateView()

    def setUnits(self):
        self.controller.updateView()

    def calcCycle(self):
        self.controller.calc()

if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    mw = MainWindow()
    mw.setWindowTitle('Air Standard Cycle Calculator')
    sys.exit(app.exec())
