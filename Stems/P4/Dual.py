# Dual.py
#region imports
from Air import *
from matplotlib import pyplot as plt
from PyQt5 import QtWidgets as qtw
import numpy as np
import sys
#endregion

class dualCycleModel:
    """
    Air standard dual (mixed) cycle:
    1->2 Isentropic compression
    2->3 Constant-volume heat addition
    3->4 Constant-pressure heat addition
    4->5 Isentropic expansion
    5->1 Constant-volume heat rejection
    """
    def __init__(self,
                 p_initial=1E5,
                 v_cylinder=3E-3,
                 t_initial=300.0,
                 ratio=18.0,
                 pressure_ratio=1.5,
                 cutoff=1.2,
                 name='Air Standard Dual Cycle'):
        self.units = units()
        self.units.SI = False
        self.air = air()
        # store parameters
        self.p_initial = p_initial
        self.T_initial = t_initial
        self.V_Cylinder = v_cylinder
        self.Ratio = ratio
        self.PressureRatio = pressure_ratio
        self.Cutoff = cutoff
        self.name = name
        # calculate number of moles in cylinder
        self.air.set(P=self.p_initial, T=self.T_initial)
        self.air.n = self.V_Cylinder / self.air.State.v
        self.air.m = self.air.n * self.air.MW
        # compute states and performance
        self._calc_states()
        self._calc_performance()
        # prepare curves for plotting
        self.upperCurve = StateDataForPlotting()
        self.lowerCurve = StateDataForPlotting()
        self.calculated = True
        self.cycleType = 'dual'

    def _calc_states(self):
        # State 1: initial
        self.State1 = self.air.set(P=self.p_initial, T=self.T_initial, name='State 1')
        # 1->2: isentropic compression
        self.State2 = self.air.set(v=self.State1.v / self.Ratio, s=self.State1.s, name='State 2')
        # 2->3: constant-volume heat addition
        P3 = self.State2.P * self.PressureRatio
        self.State3 = self.air.set(P=P3, v=self.State2.v, name='State 3')
        # 3->4: constant-pressure heat addition
        self.State4 = self.air.set(P=self.State3.P, v=self.State3.v * self.Cutoff, name='State 4')
        # 4->5: isentropic expansion
        self.State5 = self.air.set(v=self.State1.v, s=self.State4.s, name='State 5')
        # 5->1: constant-volume heat rejection
        self.State6 = self.air.set(v=self.State1.v, T=self.T_initial, name='State 1')

    def _calc_performance(self):
        n = self.air.n
        # work and heat terms
        W_comp = n * (self.State2.u - self.State1.u)
        W_exp  = n * (self.State4.u - self.State5.u)
        Q_v     = n * (self.State3.u - self.State2.u)
        Q_p     = n * (self.State4.h - self.State3.h)
        self.Q_In  = Q_v + Q_p
        self.Q_Out = n * (self.State5.u - self.State1.u)
        self.W_Cycle = W_exp - W_comp
        self.W_Compression = W_comp
        self.W_Expansion = W_exp
        self.Eff = 100.0 * self.W_Cycle / self.Q_In if self.Q_In else 0.0

    def getSI(self):
        return self.units.SI

class dualCycleController:
    def __init__(self, model=None, ax=None):
        self.model = dualCycleModel() if model is None else model
        self.view = dualCycleView()
        self.view.ax = ax

    def setWidgets(self, w):
        (self.view.lbl_THigh, self.view.lbl_TLow, self.view.lbl_P0, self.view.lbl_V0, self.view.lbl_CR,
         self.view.le_THigh, self.view.le_TLow, self.view.le_P0, self.view.le_V0, self.view.le_CR,
         self.view.lbl_PR, self.view.le_PR, self.view.lbl_CO, self.view.le_CO,
         self.view.le_T1, self.view.le_T2, self.view.le_T3, self.view.le_T4,
         self.view.lbl_T1Units, self.view.lbl_T2Units, self.view.lbl_T3Units, self.view.lbl_T4Units,
         self.view.le_PowerStroke, self.view.le_CompressionStroke, self.view.le_HeatAdded, self.view.le_Efficiency,
         self.view.lbl_PowerStrokeUnits, self.view.lbl_CompressionStrokeUnits, self.view.lbl_HeatInUnits,
         self.view.rdo_Metric, self.view.cmb_Abcissa, self.view.cmb_Ordinate,
         self.view.chk_LogAbcissa, self.view.chk_LogOrdinate, self.view.ax, self.view.canvas) = w

    def calc(self):
        T0 = float(self.view.le_TLow.text())
        P0 = float(self.view.le_P0.text())
        V0 = float(self.view.le_V0.text())
        CR = float(self.view.le_CR.text())
        PR = float(self.view.le_PR.text())
        CO = float(self.view.le_CO.text())
        SI = self.view.rdo_Metric.isChecked()
        self.set(T_0=T0, P_0=P0, V_0=V0, ratio=CR, pressure_ratio=PR, cutoff=CO, SI=SI)

    def set(self, T_0=300.0, P_0=1.0, V_0=1.0, ratio=18.0, pressure_ratio=1.5, cutoff=1.2, SI=True):
        m = self.model
        m.units.set(SI=SI)
        m.T_initial = T_0 if SI else T_0 / m.units.CF_T
        m.p_initial = P_0 if SI else P_0 / m.units.CF_P
        m.V_Cylinder = V_0 if SI else V_0 / m.units.CF_V
        m.Ratio = ratio
        m.PressureRatio = pressure_ratio
        m.Cutoff = cutoff
        m._calc_states()
        m._calc_performance()
        m.calculated = True
        self.buildDataForPlotting()
        self.updateView()

    def buildDataForPlotting(self):
        m = self.model
        a = air()
        m.lowerCurve.clear()
        for v in np.linspace(m.State1.v, m.State2.v, 30):
            s = a.set(v=v, s=m.State1.s)
            m.lowerCurve.add((s.T, s.P, s.u, s.h, s.s, s.v))
        m.upperCurve.clear()
        for T in np.linspace(m.State2.T, m.State3.T, 20):
            s = a.set(P=m.State2.P * m.PressureRatio * (T/m.State2.T), v=m.State2.v)
            m.upperCurve.add((s.T, s.P, s.u, s.h, s.s, s.v))
        for v in np.linspace(m.State3.v, m.State4.v, 20):
            s = a.set(P=m.State3.P, v=v)
            m.upperCurve.add((s.T, s.P, s.u, s.h, s.s, s.v))
        for v in np.linspace(m.State4.v, m.State5.v, 20):
            s = a.set(v=v, s=m.State4.s)
            m.upperCurve.add((s.T, s.P, s.u, s.h, s.s, s.v))
        for T in np.linspace(m.State5.T, m.State6.T, 20):
            s = a.set(v=m.State1.v, T=m.T_initial)
            m.upperCurve.add((s.T, s.P, s.u, s.h, s.s, s.v))

    def plot_cycle_XY(self, X='s', Y='T', logx=False, logy=False, mass=False, total=False):
        self.view.plot_cycle_XY(self.model, X=X, Y=Y, logx=logx, logy=logy, mass=mass, total=total)

    def updateView(self):
        self.view.updateView(self.model)

class dualCycleView:
    def __init__(self):
        # placeholders for all widgets
        self.lbl_THigh = qtw.QLabel()
        self.lbl_TLow  = qtw.QLabel()
        self.lbl_P0    = qtw.QLabel()
        self.lbl_V0    = qtw.QLabel()
        self.lbl_CR    = qtw.QLabel()
        self.lbl_PR    = qtw.QLabel()
        self.lbl_CO    = qtw.QLabel()
        self.le_THigh  = qtw.QLineEdit()
        self.le_TLow   = qtw.QLineEdit()
        self.le_P0     = qtw.QLineEdit()
        self.le_V0     = qtw.QLineEdit()
        self.le_CR     = qtw.QLineEdit()
        self.le_PR     = qtw.QLineEdit()
        self.le_CO     = qtw.QLineEdit()
        self.le_T1     = qtw.QLineEdit()
        self.le_T2     = qtw.QLineEdit()
        self.le_T3     = qtw.QLineEdit()
        self.le_T4     = qtw.QLineEdit()
        self.lbl_T1Units = qtw.QLabel()
        self.lbl_T2Units = qtw.QLabel()
        self.lbl_T3Units = qtw.QLabel()
        self.lbl_T4Units = qtw.QLabel()
        self.le_PowerStroke       = qtw.QLineEdit()
        self.le_CompressionStroke = qtw.QLineEdit()
        self.le_HeatAdded         = qtw.QLineEdit()
        self.le_Efficiency        = qtw.QLineEdit()
        self.lbl_PowerStrokeUnits = qtw.QLabel()
        self.lbl_CompressionStrokeUnits = qtw.QLabel()
        self.lbl_HeatInUnits           = qtw.QLabel()
        self.rdo_Metric  = qtw.QRadioButton()
        self.cmb_Abcissa = qtw.QComboBox()
        self.cmb_Ordinate= qtw.QComboBox()
        self.chk_LogAbcissa = qtw.QCheckBox()
        self.chk_LogOrdinate = qtw.QCheckBox()
        self.ax     = None
        self.canvas = None

    def updateView(self, cycle):
        cycle.units.set(SI=self.rdo_Metric.isChecked())
        logx = self.chk_LogAbcissa.isChecked()
        logy = self.chk_LogOrdinate.isChecked()
        xvar = self.cmb_Abcissa.currentText()
        yvar = self.cmb_Ordinate.currentText()
        if cycle.calculated:
            self.plot_cycle_XY(cycle, X=xvar, Y=yvar, logx=logx, logy=logy, mass=False, total=True)
        self.updateDisplayWidgets(cycle)

    def plot_cycle_XY(self, cycle, X='s', Y='T', logx=False, logy=False, mass=False, total=False):
        # reuse Diesel plotting
        from Diesel import dieselCycleView
        viewer = dieselCycleView()
        viewer.ax = self.ax
        viewer.convertDataCol = lambda c,data,colName,mass,total: dieselCycleView().convertDataCol(cycle,data,colName,mass,total)
        viewer.plot_cycle_XY(cycle, X, Y, logx, logy, mass, total)

    def updateDisplayWidgets(self, Model):
        U = Model.units
        SI = U.SI
        self.lbl_THigh.setText(f'T High ({U.TUnits})')
        self.lbl_TLow.setText(f'T Low ({U.TUnits})')
        self.lbl_P0.setText(f'P0 ({U.PUnits})')
        self.lbl_V0.setText(f'V0 ({U.VUnits})')
        self.lbl_CR.setText('Compression Ratio')
        self.lbl_PR.setText('Pressure Ratio')
        self.lbl_CO.setText('Cutoff Ratio')
        CFP, CFV = (1.0,1.0) if SI else (U.CF_P, U.CF_V)
        if Model.calculated:
            self.le_TLow.setText(f'{Model.T_initial if SI else U.T_KtoR(Model.T_initial):.2f}')
            self.le_P0.setText(f'{Model.p_initial*CFP:.2f}')
            self.le_V0.setText(f'{Model.V_Cylinder*CFV:.4f}')
            self.le_CR.setText(f'{Model.Ratio:.2f}')
            self.le_PR.setText(f'{Model.PressureRatio:.2f}')
            self.le_CO.setText(f'{Model.Cutoff:.2f}')
            CFE = 1.0 if SI else U.CF_E
            self.le_Efficiency.setText(f'{Model.Eff:.3f}')
            self.le_PowerStroke.setText(f'{Model.air.n*Model.W_Expansion*CFE:.3f}')
            self.le_CompressionStroke.setText(f'{Model.air.n*Model.W_Compression*CFE:.3f}')
            self.le_HeatAdded.setText(f'{Model.air.n*Model.Q_In*CFE:.3f}')
        Model.units.changed = False

# test runner
if __name__ == '__main__':
    app = qtw.QApplication(sys.argv)
    dc = dualCycleController()
    dc.set(T_0=300, P_0=1, V_0=0.003, ratio=18, pressure_ratio=1.5, cutoff=1.2, SI=False)
    dc.plot_cycle_XY(X='v', Y='P', total=True)
    sys.exit(app.exec())
