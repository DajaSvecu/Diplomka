import sys
from PyQt5.QtGui import QFont, QDoubleValidator  # type: ignore
from PyQt5.QtWidgets import (QFileDialog, QHBoxLayout, QVBoxLayout, QCheckBox,  # type: ignore
    QLabel, QLineEdit, QPushButton, QGroupBox, QRadioButton, QComboBox, QMessageBox,
    QFormLayout, QDialogButtonBox, QApplication, QDialog, QMainWindow)

from heat_exchanger_ui import Ui_MainWindow
from heat_exchanger_inputs import UserInputMedium, UserInputRest, UserInputShell
import entrance
from Medium import Medium

class HeatExchangerWindow(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.input_tube = UserInputMedium('Tube')
        #self.input_shell = UserInputMedium('Shell') ve finalni verzi
        self.input_shell = UserInputShell('Shell') # pouze docasny
        self.input_rest = UserInputRest()
        self.add_medium_inputs(self.verticalLayoutTube, self.input_tube)
        self.add_medium_inputs(self.verticalLayoutShell, self.input_shell)
        self.add_main_inputs(self.verticalLayoutRest, self.input_rest)
        self.buttonRun.pressed.connect(lambda: self.run_simulation())

    def add_medium_inputs(self, parent_layout, input) -> None:
        self.add_main_inputs(parent_layout, input)
        title_layout = QHBoxLayout()
        for text in ['Medium', 'Procenta']:
            title = QLabel(text)
            title_layout.addWidget(title)
        parent_layout.addLayout(title_layout)
        for i in range(1):
            new_layout = QHBoxLayout()
            for parameter in input.medium:
                name = input.name + parameter['name'] + str(i)
                setattr(self, name, QLineEdit())
                getattr(self, name).setText(str(parameter['default_value']))
                #getattr(self, name).setDisabled(True)
                new_layout.addWidget(getattr(self, name))
            parent_layout.addLayout(new_layout)

    def add_main_inputs(self, parent_layout, input) -> None:
        title_text = '{}'.format(input.name)
        title = QLabel(title_text)
        parent_layout.addWidget(title)
        
        for element in input.parameters:
            new_layout = QHBoxLayout()

            label_text = "{} [{}]:".format(element['name'], element['unit'])
            label = QLabel(label_text)
            name = input.name + element['name']
            setattr(self, name, QLineEdit())
            getattr(self, name).setText(str(element['default_value']))
            new_layout.addWidget(label)
            new_layout.addWidget(getattr(self, name))

            parent_layout.addLayout(new_layout)
        
    def get_medium_inputs(self, input) -> dict:
        values = self.get_main_inputs(input)
        medium = [[0,0] for y in range(1)]
        for i in range(1):
            j = 0
            for element in input.medium:
                name = input.name + element['name'] + str(i)
                value = getattr(self, name).text()
                parse_value = element['parse_function'](value)
                medium[i][j] = parse_value
                j += 1
        values['Medium'] = medium
        return values

    def get_main_inputs(self, input) -> dict:
        values = {}
        for element in input.parameters:
            name = input.name + element['name']
            value = getattr(self, name).text()
            parse_value = element['parse_function'](value)
            values[element['name']] = parse_value
        return values
    
    
    def run_simulation(self):
        print('RUNNING')
        medium_tube = self.get_medium_inputs(self.input_tube)
        medium_shell = self.get_medium_inputs(self.input_shell)
        rest = self.get_main_inputs(self.input_rest)
        print(rest)
        #result = calculate_everything(medium_tube, medium_shell, rest)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    heat_exchanger = HeatExchangerWindow()
    heat_exchanger.show()
    sys.exit(app.exec_())