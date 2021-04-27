# Created by Hao JIN on 2021/04/27.

import sys
import os
import json
import numpy as np

from PyQt6.QtWidgets import QApplication, QMainWindow, QLabel
from PyQt6.QtGui import QIcon, QPixmap
from impulse_ui import Ui_MainWindow


class MyMainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super().__init__()

        # set up the user interface from Designer
        self.setupUi(self)

        # make some local modifications
        self.readJson()

        # connect up the buttons
        self.pushButton.clicked.connect(self.design)
        self.comboBox.currentTextChanged.connect(self.itemSelected)
    
    # comboBox
    def itemSelected(self):
        item = self.comboBox.currentText()
        self.lineEdit_v.setText(str(self.materials[item]["velocity"]))
        self.lineEdit_K2.setText(str(self.materials[item]["K2"]))
        self.lineEdit_Co.setText(str(self.materials[item]["Co"]))

    # read .json file
    def readJson(self):
        self.json_file = "impulse_materials.json"
        with open(self.json_file) as f_obj:
            self.materials = json.load(f_obj)
    
    # calculate the lambda
    def cal_lambda(self):
        self.lambda0 = np.round(self.v/self.freq, 1)
        self.lineEdit_lambda.setText(str(self.lambda0))
    
    # calculate the number of finger pairs
    def cal_Np(self):
        self.Np = int(2*self.freq/self.nbw)
        self.lineEdit_Np.setText(str(self.Np))
    
    def cal_aperture(self):
        temp_a = (1/self.res)*(1/(2*self.freq*1e-4*self.Co*self.Np))*(4*self.K2*self.Np)
        temp_b = ((4*self.K2*self.Np)**2 + np.pi**2)
        self.aperture = np.round(temp_a/temp_b*1e6, 1)
        self.lineEdit_aperture.setText(str(self.aperture))

    # design
    def design(self):
        self.freq = float(self.lineEdit_freq.text())
        self.nbw = float(self.lineEdit_nbw.text())
        self.res = float(self.lineEdit_res.text())
        self.v = float(self.lineEdit_v.text())
        self.K2 = float(self.lineEdit_K2.text())*0.01
        self.Co = float(self.lineEdit_Co.text())
        self.cal_lambda()
        self.cal_Np()
        self.cal_aperture()
        

app = QApplication(sys.argv)
window = MyMainWindow()
window.show()
sys.exit(app.exec())