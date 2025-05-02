# main.py
import sys
from PyQt5 import QtWidgets
from sto_view import STOView
from sto_controller import STOController

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    view = STOView()
    controller = STOController(view)
    view.show()
    sys.exit(app.exec_())