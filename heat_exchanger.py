import sys
from PyQt5.QtGui import QFont, QDoubleValidator  # type: ignore
from PyQt5.QtWidgets import (QFileDialog, QHBoxLayout, QVBoxLayout, QCheckBox,  # type: ignore
    QLabel, QLineEdit, QPushButton, QGroupBox, QRadioButton, QComboBox, QMessageBox,
    QFormLayout, QDialogButtonBox, QApplication, QDialog, QMainWindow, QTableWidgetItem)
from PyQt5.QtCore import Qt
import matplotlib.pyplot as plt

from heat_exchanger_ui import Ui_MainWindow
from heat_exchanger_inputs import UserInputMedium, UserInputRest, UserInputShell
from Medium import Medium
from heat_exchanger_calc import Calculate

class HeatExchangerWindow(QMainWindow, Ui_MainWindow):
    """
    Class obsahující komunikaci mezi uživatelským prostředím a backendem.
    Samotné uživatelské prostředí se nachází v souboru heat_exchanger_ui odkud je importováno.
    """
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
        self.prepare_table()

        self.buttonRun.pressed.connect(lambda: self.run_simulation())

    def add_medium_inputs(self, parent_layout, input) -> None:
        """
        Vygenerovani vstupu pro jednotliva media do uzivatelskeho prostredi.
        """
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
        """
        Pridani vstupu hlavni vstupu do uzivatelskeho prostredi.
        """
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
            if element['name'] == 'T2': getattr(self, name).setDisabled(True)
            new_layout.addWidget(label)
            new_layout.addWidget(getattr(self, name))

            parent_layout.addLayout(new_layout)

    def get_medium_inputs(self, input) -> dict:
        """
        Ziskani vstupu od uzivatele o mediu.
        """
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
        """
        Ziskani hlavnich vstupu od uzivatele
        """
        values = {}
        for element in input.parameters:
            name = input.name + element['name']
            value = getattr(self, name).text()
            parse_value = element['parse_function'](value)
            values[element['name']] = parse_value
        return values
    
    
    def run_simulation(self):
        """
        Spusteni po stisknuti tlacitka
        Slouzi k vypoctu bilancni rovnice (vypoctu vystupnich teplot)
        Slouzi k vypoctu jednotlivym vymeniku
        """
        print('RUNNING')
        self.table.clearContents()
        self.table.setRowCount(0)
        medium_tube = self.get_medium_inputs(self.input_tube)
        medium_shell = self.get_medium_inputs(self.input_shell)
        rest = self.get_main_inputs(self.input_rest)
        try:
            calculate = Calculate(medium_tube, medium_shell, rest)
            getattr(self, 'TubeT2').setText(str(round(calculate.tube.t2, 2)))
            getattr(self, 'ShellT2').setText(str(round(calculate.shell.t2, 2)))
        except ValueError as error:
            print(error.args)
        except Exception as error:
            print(error.args)
        try:
            vysledky = calculate.calculate_all()
            
        except:
            pass
        else:
            print('Pozadavky na vymenik splnilo {} vymeniku.'.format(len(vysledky)))
            self.show_output(vysledky)
            print('DONE!')
            for vysledek in vysledky:
                plt.scatter(vysledek[1]['hmotnost'],vysledek[1]['tlak_ztraty'])
            plt.show()    

    def prepare_table(self):
        """
        Pripraveni sloupcu v tabulce
        """
        i = 0
        for item in ['DN', 'd_in', 'tl', 't_p', 't_t']:
            self.table.insertColumn(i)
            self.table.setHorizontalHeaderItem(i, QTableWidgetItem(item))
            i += 1
        for item in ['vyska_prep', 'tl_prep', 'delka','pocet_prepazek', 'pocet_trubek', 'w_tube', 'w_shell', 'tlak_ztraty', 'hmotnost']:
            self.table.insertColumn(i)
            self.table.setHorizontalHeaderItem(i, QTableWidgetItem(item))
            i += 1
        

    def show_output(self, outputs):
        """
        Vygenerovani vymeniku ktere splnili pozadavky
        """
        i = 0
        self.table.setSortingEnabled(True)
        for output in outputs:
            self.table.insertRow(i)
            j = 0
            for part_output in output:
                for x in part_output:
                    item = QTableWidgetItem()
                    if x == 'shell':
                        item.setData(0, part_output[x]['DN'])
                    else:
                        item.setData(0, part_output[x])
                    item.setFlags(Qt.ItemFlags(1))
                    self.table.setItem(i, j, item)
                    j += 1
            i += 1


            
        

if __name__ == '__main__':
    app = QApplication(sys.argv)
    heat_exchanger = HeatExchangerWindow()
    heat_exchanger.show()
    sys.exit(app.exec_())