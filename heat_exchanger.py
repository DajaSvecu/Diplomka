import sys
from PyQt5.QtGui import QFont, QDoubleValidator  # type: ignore
from PyQt5.QtCore import QThreadPool, Qt  # type: ignore
from PyQt5.QtWidgets import (QFileDialog, QHBoxLayout, QVBoxLayout, QCheckBox,  # type: ignore
    QLabel, QLineEdit, QPushButton, QGroupBox, QRadioButton, QComboBox, QMessageBox,
    QFormLayout, QDialogButtonBox, QApplication, QDialog, QMainWindow)

from heat_exchanger_ui import Ui_MainWindow
from heat_exchanger_inputs import UserInput

class HeatExchangerWindow(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.input_entrance = UserInput()
        self.add_user_inputs(self.verticalLayoutEntrance)
        self.buttonRun.pressed.connect(lambda: self.run_simulation())
    def add_user_inputs(self, parent_layout) -> None:
        
        for entry in self.input_entrance.parameters_from_user:
            new_layout = QHBoxLayout()

            label_text = "{} [{}]:".format(entry['variable_name'], entry['unit'])
            label = QLabel(label_text)
            label.setToolTip(entry['description'])
            setattr(self, entry['name'], QLineEdit())
            new_layout.addWidget(label)
            new_layout.addWidget(getattr(self, entry['name']))

            parent_layout.addLayout(new_layout)
    def get_user_inputs(self) -> dict:
        values = {}
        for element in self.input_entrance.parameters_from_user:
            value = getattr(self, element['name']).text()
            parse_value = element['parse_function'](value)
            values[element['variable_name']] = parse_value
        return values
    def run_simulation(self):
        print('RUNNING')
        inputs = self.get_user_inputs()
        print(inputs)
    def entrance(self):
        pass
    def medium(self):
        pass
    def validate_input(self):
        pass

if __name__ == '__main__':
    app = QApplication(sys.argv)
    heat_exchanger = HeatExchangerWindow()
    heat_exchanger.show()
    sys.exit(app.exec_())