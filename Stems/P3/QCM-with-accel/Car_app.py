#region imports
from Car_GUI import Ui_Form
import sys
from PyQt5 import QtCore as qtc
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from QuarterCarModel import CarController
#endregion

class MainWindow(qtw.QWidget, Ui_Form):
    def __init__(self):
        """
        Main window constructor.
        """
        super().__init__()
        # call setupUi from Ui_Form
        self.setupUi(self)

        # setup car controller
        input_widgets = (
            self.le_m1, self.le_v, self.le_k1, self.le_c1,
            self.le_m2, self.le_k2, self.le_ang, self.le_tmax,
            self.chk_IncludeAccel
        )
        display_widgets = (
            self.gv_Schematic,
            self.chk_LogX, self.chk_LogY,
            self.chk_LogAccel, self.chk_ShowAccel,
            self.lbl_MaxMinInfo, self.layout_Plot
        )

        # instantiate controller
        self.controller = CarController((input_widgets, display_widgets))

        # connect signals
        self.btn_calculate.clicked.connect(self.controller.calculate)
        self.pb_Optimize.clicked.connect(self.doOptimize)
        self.chk_LogX.stateChanged.connect(self.controller.doPlot)
        self.chk_LogY.stateChanged.connect(self.controller.doPlot)
        self.chk_LogAccel.stateChanged.connect(self.controller.doPlot)
        self.chk_ShowAccel.stateChanged.connect(self.controller.doPlot)

        # install filters & canvas listeners
        self.controller.setupEventFilter(self)
        self.controller.setupCanvasMoveEvent(self)

        self.show()

    def eventFilter(self, obj, event):
        """
        Override to track mouse / wheel on the schematic.
        """
        if obj == self.gv_Schematic.scene():
            et = event.type()
            if et == qtc.QEvent.GraphicsSceneMouseMove:
                pos = event.scenePos()
                title = f"Mouse Position:  x={pos.x():.2f}, y={-pos.y():.2f}"
                self.setWindowTitle(title)
            elif et == qtc.QEvent.GraphicsSceneWheel:
                # zoom with wheel — only here do we redraw
                zm = self.controller.getZoom()
                zm += 0.1 if event.delta() > 0 else -0.1
                zm = max(0.1, zm)
                self.controller.setZoom(zm)
                self.controller.updateSchematic()
        return super(MainWindow, self).eventFilter(obj, event)

    def mouseMoveEvent_Canvas(self, event):
        """
        Display live interpolation values from the plot canvas.
        """
        if event.inaxes:
            t = event.xdata
            self.controller.animate(t)
            yw, yb, yr, acc = self.controller.getPoints(t)
            self.setWindowTitle(
                f"t={t:.2f}s, road={yr*1000:.2f}mm, wheel={yw*1000:.2f}mm, "
                f"body={yb*1000:.2f}mm, accel={acc:.2f}g"
            )

    def doOptimize(self):
        """
        Run the suspension optimizer.
        """
        app.setOverrideCursor(qtc.Qt.WaitCursor)
        self.controller.OptimizeSuspension()
        app.restoreOverrideCursor()


if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    mw  = MainWindow()
    mw.setWindowTitle('Quarter Car Model')
    sys.exit(app.exec_())
