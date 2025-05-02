#region imports
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import numpy as np
import math
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtCore as qtc
from PyQt5 import QtGui as qtg

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
#endregion

class CarModel():
    def __init__(self):
        # time array (updated in calculate)
        self.timeData = np.linspace(0,3,200)
        # data containers
        self.roadPosData     = []
        self.wheelPosData    = []
        self.bodyPosData     = []
        self.accelBodyData   = []
        self.springForceData = []
        self.dashpotForceData= []

        # default inputs
        self.tmax    = 3.0
        self.tramp   = 1.0
        self.yangdeg = 45.0
        self.ymag    = 6.0/(12*3.3)

        self.m1 = 450
        self.m2 = 20
        self.k1 = 15000
        self.k2 = 90000
        self.c1 = 4500
        self.v  = 120.0

        # bounds for optimizer/info
        self.mink1    = (self.m1*9.81)/(6.0*25.4/1000.0)
        self.maxk1    = (self.m1*9.81)/(3.0*25.4/1000.0)
        self.mink2    = ((self.m1+self.m2)*9.81)/(1.5*25.4/1000.0)
        self.maxk2    = ((self.m1+self.m2)*9.81)/(0.75*25.4/1000.0)
        self.accelLim = 1.5
        self.SSE      = 0.0

class CarView():
    def __init__(self, model, args):
        self.model = model
        # unpack widgets
        self.input_widgets, self.display_widgets = args
        ( self.le_m1, self.le_v, self.le_k1, self.le_c1,
          self.le_m2, self.le_k2, self.le_ang, self.le_tmax,
          self.chk_IncludeAccel ) = self.input_widgets

        ( self.gv_Schematic,
          self.chk_LogX, self.chk_LogY,
          self.chk_LogAccel, self.chk_ShowAccel,
          self.lbl_MaxMinInfo, self.layout_Plot ) = self.display_widgets

        # create the tab widget
        self.tabs = qtw.QTabWidget()
        self.layout_Plot.addWidget(self.tabs)

        # — Position vs Time Tab —
        self.tab_pos = qtw.QWidget()
        self.tabs.addTab(self.tab_pos, "Position vs. Time")
        pos_layout = qtw.QVBoxLayout(self.tab_pos)
        self.fig_pos    = Figure(tight_layout=True)
        self.canvas_pos = FigureCanvasQTAgg(self.fig_pos)
        pos_layout.addWidget(NavigationToolbar2QT(self.canvas_pos))
        pos_layout.addWidget(self.canvas_pos)
        self.ax_pos = self.fig_pos.add_subplot(111)
        self.ax_acc = self.ax_pos.twinx()

        # — Force vs Time Tab —
        self.tab_force = qtw.QWidget()
        self.tabs.addTab(self.tab_force, "Force vs. Time")
        f_layout = qtw.QVBoxLayout(self.tab_force)
        self.fig_force    = Figure(tight_layout=True)
        self.canvas_force = FigureCanvasQTAgg(self.fig_force)
        f_layout.addWidget(NavigationToolbar2QT(self.canvas_force))
        f_layout.addWidget(self.canvas_force)
        self.ax_force = self.fig_force.add_subplot(111)

        # build initial schematic
        self.buildScene()

    def buildScene(self):
        self.scene = qtw.QGraphicsScene()
        self.scene.setSceneRect(-200, -200, 400, 400)
        self.gv_Schematic.setScene(self.scene)
        self.setupPensAndBrushes()
        # … instantiate your Mass, Wheel, Spring, Dashpot, Road items here …

    def setupPensAndBrushes(self):
        self.penTire  = qtg.QPen(qtc.Qt.black); self.penTire.setWidth(3)
        self.penMass  = qtg.QPen(qtc.Qt.black); self.penMass.setWidth(1)
        color = qtg.QColor(qtc.Qt.gray); color.setAlpha(64)
        self.brushWheel= qtg.QBrush(color)
        self.brushMass = qtg.QBrush(qtg.QColor(200,200,200,64))
        self.groundPen   = qtg.QPen(qtc.Qt.black); self.groundPen.setWidth(1)
        self.groundBrush = qtg.QBrush(qtc.Qt.black); self.groundBrush.setStyle(qtc.Qt.DiagCrossPattern)

    def setupEventFilter(self, window):
        self.gv_Schematic.setMouseTracking(True)
        self.gv_Schematic.scene().installEventFilter(window)

    def setupCanvasMoveEvent(self, window):
        self.canvas_pos.mpl_connect("motion_notify_event", window.mouseMoveEvent_Canvas)

    def updateView(self, model):
        # push model → UI
        self.le_m1.setText(f"{model.m1:.2f}")
        self.le_k1.setText(f"{model.k1:.2f}")
        self.le_c1.setText(f"{model.c1:.2f}")
        self.le_m2.setText(f"{model.m2:.2f}")
        self.le_k2.setText(f"{model.k2:.2f}")
        self.le_ang.setText(f"{model.yangdeg:.2f}")
        self.le_tmax.setText(f"{model.tmax:.2f}")
        info = (
            f"k1_min={model.mink1:.2f}, k1_max={model.maxk1:.2f}\n"
            f"k2_min={model.mink2:.2f}, k2_max={model.maxk2:.2f}\n"
            f"SSE={model.SSE:.2f}"
        )
        self.lbl_MaxMinInfo.setText(info)
        self.doPlot(model)

    def doPlot(self, model):
        t = model.timeData

        # — POSITION PLOT —
        self.ax_pos.clear(); self.ax_acc.clear()
        if self.chk_LogX.isChecked():     self.ax_pos.set_xscale('log')
        if self.chk_LogY.isChecked():     self.ax_pos.set_yscale('log')
        if self.chk_LogAccel.isChecked(): self.ax_acc.set_yscale('log')

        self.ax_pos.plot(t, model.bodyPosData,  label='Body Position')
        self.ax_pos.plot(t, model.wheelPosData, label='Wheel Position')
        self.ax_pos.plot(t, model.roadPosData,  color='k', lw=3, label='Road')
        self.ax_pos.set_xlabel('time (s)')
        self.ax_pos.set_ylabel('Position (m)')
        self.ax_pos.legend(loc='upper left')

        if self.chk_ShowAccel.isChecked():
            self.ax_acc.plot(t, model.accelBodyData, 'g-', label='Acceleration (g)')
            self.ax_acc.legend(loc='upper right')

        self.canvas_pos.draw()

        # — FORCE PLOT —
        self.ax_force.clear()
        self.ax_force.plot(t, model.springForceData,  label='Spring Force (N)')
        self.ax_force.plot(t, model.dashpotForceData, label='Dashpot Force (N)')
        self.ax_force.set_xlabel('time (s)')
        self.ax_force.set_ylabel('Force (N)')
        self.ax_force.legend()
        self.canvas_force.draw()

class CarController():
    def __init__(self, args):
        self.input_widgets, self.display_widgets = args
        self.model = CarModel()
        self.view  = CarView(self.model, args)

    def setupEventFilter(self, window):
        self.view.setupEventFilter(window)

    def setupCanvasMoveEvent(self, window):
        self.view.setupCanvasMoveEvent(window)

    def ode_system(self, t, X):
        x1, x1dot, x2, x2dot = X
        # ramp‐to‐ymag road input
        if t < self.model.tramp:
            y = self.model.ymag * (t/self.model.tramp)
        else:
            y = self.model.ymag

        x1dd = (1/self.model.m1)*( self.model.c1*(x2dot - x1dot)
                                  + self.model.k1*(x2 - x1) )
        x2dd = (1/self.model.m2)*( -self.model.c1*(x2dot - x1dot)
                                  - self.model.k1*(x2 - x1)
                                  + self.model.k2*(y - x2) )
        return [x1dot, x1dd, x2dot, x2dd]

    def calculate(self, doCalc=True):
        # read UI → model
        self.model.m1      = float(self.view.le_m1.text())
        self.model.m2      = float(self.view.le_m2.text())
        self.model.c1      = float(self.view.le_c1.text())
        self.model.k1      = float(self.view.le_k1.text())
        self.model.k2      = float(self.view.le_k2.text())
        self.model.v       = float(self.view.le_v.text())
        self.model.yangdeg = float(self.view.le_ang.text())
        self.model.tmax    = float(self.view.le_tmax.text())

        # derived parameters
        v_ms = 1000 * self.model.v / 3600
        self.model.angrad  = math.radians(self.model.yangdeg)
        self.model.ymag    = 6.0/(12*3.3)
        self.model.tramp   = self.model.ymag/(math.sin(self.model.angrad)*v_ms)
        self.model.timeData = np.logspace(-6, np.log10(self.model.tmax), 2000)

        # solve ODE
        sol = solve_ivp(
            self.ode_system,
            [0, self.model.tmax],
            [0,0,0,0],
            t_eval=self.model.timeData
        )

        # unpack
        self.model.roadPosData  = [
            (self.model.ymag*(t/self.model.tramp) if t < self.model.tramp else self.model.ymag)
            for t in sol.t
        ]
        self.model.bodyPosData  = sol.y[0]
        self.model.wheelPosData = sol.y[2]

        # acceleration (g)
        vel   = sol.y[1]
        dt    = np.diff(sol.t)
        accel = np.empty_like(sol.t)
        accel[:-1] = np.diff(vel)/(9.81*dt)
        accel[-1]  = accel[-2]
        self.model.accelBodyData = accel

        # ** compute forces **
        self.model.springForceData  = self.model.k1*(self.model.wheelPosData - self.model.bodyPosData)
        self.model.dashpotForceData = self.model.c1*(sol.y[3] - sol.y[1])

        # update UI
        self.view.updateView(self.model)

    def doPlot(self):
        # re-plot data without re-solving
        self.view.doPlot(self.model)

    # stubs for schematic & animation
    def getZoom(self):
        return self.view.gv_Schematic.transform().m11()

    def setZoom(self, z):
        self.view.gv_Schematic.resetTransform()
        self.view.gv_Schematic.scale(z, z)

    def updateSchematic(self):
        self.view.buildScene()

    def animate(self, t):
        pass

    def getPoints(self, t):
        return (0, 0, 0, 0)

    def OptimizeSuspension(self):
        pass
