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
        
        self.input_tube = UserInputMedium('Tube', 'Trubkovy prostor')
        self.input_shell = UserInputShell('Shell', 'Mezitrubkovy prostor') # pouze docasny na testovani
        #self.input_shell = UserInputMedium('Shell', 'Mezitrubkovy prostor') #ve finalni verzi
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
        for parameter in input.medium:
            title = QLabel(parameter['name'])
            title.setFont(QFont("Times", 12, QFont.Bold))
            title.setToolTip(parameter['hint'])
            title_layout.addWidget(title)
        parent_layout.addLayout(title_layout)

        for i in range(1):
            new_layout = QHBoxLayout()
            for parameter in input.medium:
                name = input.name + parameter['name'] + str(i)
                setattr(self, name, QLineEdit())
                getattr(self, name).setFont(QFont("Times", 12))
                getattr(self, name).setText(str(parameter['default_value']))
                #getattr(self, name).setDisabled(True)
                new_layout.addWidget(getattr(self, name))
            parent_layout.addLayout(new_layout)

    def add_main_inputs(self, parent_layout, input) -> None:
        """
        Pridani vstupu hlavni vstupu do uzivatelskeho prostredi.
        """
        title = QLabel(input.title)
        title.setAlignment(Qt.AlignCenter)
        title.setFont(QFont("Times", 16, QFont.Bold))
        parent_layout.addWidget(title)
        
        for element in input.parameters:
            new_layout = QHBoxLayout()

            label_text = "{} [{}]:".format(element['name'], element['unit'])
            label = QLabel(label_text)
            label.setFont(QFont("Times", 12))
            label.setToolTip(element['hint'])
            name = input.name + element['name']
            setattr(self, name, QLineEdit())
            getattr(self, name).setText(str(element['default_value']))
            getattr(self, name).setFont(QFont("Times", 12))
            if element['name'] == 'T2': getattr(self, name).setDisabled(True)
            new_layout.addWidget(label)
            new_layout.addWidget(getattr(self, name))

            parent_layout.addLayout(new_layout)

    def get_medium_inputs(self, input) -> dict:
        """
        Ziskani vstupu od uzivatele o mediu.
        """
        values = self.get_main_inputs(input)
        medium = []
        for i in range(1):
            part_medium = [0, 0]
            j = 0
            for element in input.medium:

                name = input.name + element['name'] + str(i)
                value = getattr(self, name).text().replace(',', '.')
                if value == '':
                    break
                try:
                    parse_value = element['parse_function'](value) * element['to_SI']
                except ValueError:
                    self.show_error_dialog_to_user('Hodnota {} ma spatny format. Pro "medium" zadejte retezec a pro "procenta" cislo!'.format(name))
                
                part_medium[j] = parse_value
                j += 1
            if value != '': medium.append(part_medium)
        values['Medium'] = medium
        return values

    def get_main_inputs(self, input) -> dict:
        """
        Ziskani hlavnich vstupu od uzivatele
        """
        values = {}
        for element in input.parameters:
            try:
                name = input.name + element['name']
                value = getattr(self, name).text().replace(",", ".")
                parse_value = element['parse_function'](value) * element['to_SI']
                if parse_value < 0: raise Exception('Hodnota {} je mensi nez nula. Zadejte kladne cislo.'.format(name))
                values[element['name']] = parse_value
            except ValueError:
                self.show_error_dialog_to_user('Nezadali jste cislo u hodnoty {}!'.format(name))
            except Exception as error:
                self.show_error_dialog_to_user(error.args[0])

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
            vysledky = calculate.calculate_all()   
        except Exception as error:
            self.show_error_dialog_to_user(error.args[0])
        else:
            print('Pozadavky na vymenik splnilo {} vymeniku.'.format(len(vysledky)))
            self.show_output(vysledky)
            for vysledek in vysledky:
                plt.scatter(vysledek[1]['hmotnost'],vysledek[1]['tlak_ztraty'])
            plt.title('Graf doporucenych vymeniku')
            plt.xlabel('Hmotnost [kg]')
            plt.ylabel('Tlakove ztraty [Pa]')
            plt.show()    
            print('DONE!')

    def show_error_dialog_to_user(self, error_message: str) -> None:
        """
        Displays a separate dialog (alert) informing user of something bad, like
            invalid user input or simulation errors.

        Args:
            error_message ... what should be shown to the user
        """

        print(error_message)

        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setText("Error")
        msg.setWindowTitle("Error")
        msg.setInformativeText(error_message)
        msg.exec_()

    def prepare_table(self):
        """
        Pripraveni sloupcu v tabulce
        """
        i = 0
        for item in ['DN[-]', 'd_out[mm]', 'tl_trub[mm]', 'roztec_trub[mm]', 'delka[mm]', 'roztec_prep[mm]', 'vyska_prep[mm]']:
            self.table.insertColumn(i)
            self.table.setHorizontalHeaderItem(i, QTableWidgetItem(item))
            i += 1
        for item in ['tl_prep[mm]','pocet_prep[-]', 'pocet_trub[-]', 'TP[m/s]', 'MZP[m/s]', 'vykon [W]',
        'tlak_ztraty[Pa]', 'hmotnost[kg]']:
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
            for x in output[0]:
                item = QTableWidgetItem()
                if x == 'shell':
                    item.setData(0, output[0][x]['DN'])
                else:
                    item.setData(0, output[0][x]*1000)
                item.setFlags(Qt.ItemFlags(1))
                self.table.setItem(i, j, item)
                j += 1
            for y in output[1]:
                item = QTableWidgetItem()
                if y == 'tl_prep':
                    item.setData(0, output[1][y]*1000)
                else:
                    item.setData(0, output[1][y])
                item.setFlags(Qt.ItemFlags(1))
                self.table.setItem(i, j, item)
                j += 1
            i += 1


            
        

if __name__ == '__main__':
    app = QApplication(sys.argv)
    heat_exchanger = HeatExchangerWindow()
    heat_exchanger.show()
    sys.exit(app.exec_())