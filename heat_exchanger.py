import sys
from PyQt5.QtGui import QFont, QDoubleValidator  # type: ignore
from PyQt5.QtCore import QThreadPool, Qt  # type: ignore
from PyQt5.QtWidgets import (QFileDialog, QHBoxLayout, QVBoxLayout, QCheckBox,  # type: ignore
    QLabel, QLineEdit, QPushButton, QGroupBox, QRadioButton, QComboBox, QMessageBox,
    QFormLayout, QDialogButtonBox, QApplication, QDialog, QMainWindow)

from heat_exchanger_ui import Ui_MainWindow
from heat_exchanger_inputs import UserInputMedium, UserInputRest
import entrance
import Medium

class HeatExchangerWindow(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.input_tube = UserInputMedium('Tube')
        self.input_shell = UserInputMedium('Shell')
        self.input_rest = UserInputRest()
        self.add_user_inputs(self.verticalLayoutTube, self.input_tube)
        self.add_user_inputs(self.verticalLayoutShell, self.input_shell)
        self.add_user_inputs(self.verticalLayoutRest, self.input_rest)
        self.buttonRun.pressed.connect(lambda: self.run_simulation())


    def add_user_inputs(self, parent_layout, input) -> None:
        title_text = '{}'.format(input.name)
        title = QLabel(title_text)
        parent_layout.addWidget(title)

        for element in input.parameters:
            new_layout = QHBoxLayout()

            label_text = "{} [{}]:".format(element['name'], element['unit'])
            label = QLabel(label_text)
            name = input.name + element['name']
            setattr(self, name, QLineEdit())
            new_layout.addWidget(label)
            new_layout.addWidget(getattr(self, name))

            parent_layout.addLayout(new_layout)


    def get_user_inputs(self, input) -> dict:
        values = {}
        for element in input.parameters:
            name = input.name + element['name']
            value = getattr(self, name).text()
            parse_value = element['parse_function'](value)
            values[element['name']] = parse_value
        return values
    
    def run_simulation(self):
        print('RUNNING')
        medium_tube = self.get_user_inputs(self.input_tube)
        medium_shell = self.get_user_inputs(self.input_shell)
        Q = entrance.finish_inputs(medium_tube, medium_shell)
        Medium.props(medium_tube)
        Medium.props(medium_shell)
        print(Q)
        print(medium_shell)
        print(medium_tube)
        #self.rest = self.get_user_inputs(self.input_rest)
        

if __name__ == '__main__':
    app = QApplication(sys.argv)
    heat_exchanger = HeatExchangerWindow()
    heat_exchanger.show()
    sys.exit(app.exec_())