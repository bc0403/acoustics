# Created by Hao JIN on 2021/04/27.

import sys
import os
import numpy as np
import json
from sympy import Matrix, cos, sin, Rational, pprint, ones

from PyQt5.QtWidgets import (QApplication, QMainWindow, QLabel,
    QTableWidgetItem)
from PyQt5.QtGui import QIcon, QPixmap

from euler_ui import Ui_MainWindow

from acoustics import LN_auld as LN
from acoustics import PiezoMaterial, Trig3m, Trig32, Hex6mm
from acoustics import Trig3m, Trig32

class MyMainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super().__init__()

        # set up the user interface from Designer
        self.setupUi(self)

        # make some local modifications
        self.readJson()

        for i in range(6):
            for j in range(6):
                self.table_cE.setItem(i, j, QTableWidgetItem(str("0.0")))
                self.table_cE2.setItem(i, j, QTableWidgetItem(str("0.0")))
        
        for i in range(3):
            for j in range(6):
                self.table_e.setItem(i, j, QTableWidgetItem(str("0.0")))
                self.table_e2.setItem(i, j, QTableWidgetItem(str("0.0")))
        
        for i in range(3):
            for j in range(3):
                self.table_eS.setItem(i, j, QTableWidgetItem(str("0.0")))
                self.table_eS2.setItem(i, j, QTableWidgetItem(str("0.0")))


        # connect up the buttons
        self.pushButton.clicked.connect(self.rotate)
        self.comboBox.currentTextChanged.connect(self.itemSelected)
    
    # read .json file
    def readJson(self):
        self.json_file = "euler_materials.json"
        with open(self.json_file) as f_obj:
            self.materials = json.load(f_obj)
    
    def itemSelected(self):
        item = self.comboBox.currentText()
        if item == "LiNbO_3":  # Trig3m
            self.material = Trig3m(self.materials[item]["c11"],
                                   self.materials[item]["c12"],
                                   self.materials[item]["c13"],
                                   self.materials[item]["c14"],
                                   self.materials[item]["c33"],
                                   self.materials[item]["c44"],
                                   self.materials[item]["ex5"],
                                   self.materials[item]["ey2"],
                                   self.materials[item]["ez1"],
                                   self.materials[item]["ez3"],
                                   self.materials[item]["eSxx"],
                                   self.materials[item]["eSzz"],
                                   )
            for i in range(6):
                for j in range(6):
                    self.table_cE.setItem(i, j, QTableWidgetItem(str(float(self.material.c[i, j]))))
            
            for i in range(3):
                for j in range(6):
                    self.table_e.setItem(i, j, QTableWidgetItem(str(float(self.material.e[i, j]))))
            
            for i in range(3):
                for j in range(3):
                    self.table_eS.setItem(i, j, QTableWidgetItem(str(float(self.material.eS[i, j]))))
        
        elif item == "Quartz":  # Trig32
            self.material = Trig32(self.materials[item]["c11"],
                                   self.materials[item]["c12"],
                                   self.materials[item]["c13"],
                                   self.materials[item]["c14"],
                                   self.materials[item]["c33"],
                                   self.materials[item]["c44"],
                                   self.materials[item]["ex1"],
                                   self.materials[item]["ex4"],
                                   self.materials[item]["eSxx"],
                                   self.materials[item]["eSzz"],
                                   )
            for i in range(6):
                for j in range(6):
                    self.table_cE.setItem(i, j, QTableWidgetItem(str(float(self.material.c[i, j]))))
            
            for i in range(3):
                for j in range(6):
                    self.table_e.setItem(i, j, QTableWidgetItem(str(float(self.material.e[i, j]))))
            
            for i in range(3):
                for j in range(3):
                    self.table_eS.setItem(i, j, QTableWidgetItem(str(float(self.material.eS[i, j]))))
        
        else:
            for i in range(6):
                for j in range(6):
                    self.table_cE.setItem(i, j, QTableWidgetItem(str("0.0")))
            
            for i in range(3):
                for j in range(6):
                    self.table_e.setItem(i, j, QTableWidgetItem(str("0.0")))
            
            for i in range(3):
                for j in range(3):
                    self.table_eS.setItem(i, j, QTableWidgetItem(str("0.0")))
    
    def rotate(self):
        self.c = ones(6, 6)
        self.e = ones(3, 6)
        self.eS = ones(3, 3)

        for i in range(6):
            for j in range(6):
                self.c[i,j] = float(self.table_cE.item(i, j).text())

        for i in range(3):
            for j in range(6):
                self.e[i,j] = float(self.table_e.item(i, j).text())
        
        for i in range(3):
            for j in range(3):
                self.eS[i,j] = float(self.table_eS.item(i, j).text())
        
        self.alpha = float(self.lineEdit_alpha.text())/180*np.pi
        self.beta = float(self.lineEdit_beta.text())/180*np.pi
        self.gamma = float(self.lineEdit_gamma.text())/180*np.pi
        self.rho = 4700e3  # dummy

        self.piezoMaterial = PiezoMaterial(self.rho, self.c, self.eS, self.e)
        self.piezoMaterial.rot_euler_update(self.alpha, self.beta, self.gamma)

        for i in range(6):
            for j in range(6):
                self.table_cE2.setItem(i, j, QTableWidgetItem(str(format(float(self.piezoMaterial.stiffness[i, j]),'.1f'))))
        
        for i in range(3):
            for j in range(6):
                self.table_e2.setItem(i, j, QTableWidgetItem(str(format(float(self.piezoMaterial.piezoelec[i, j]),'.1f'))))
        
        for i in range(3):
            for j in range(3):
                self.table_eS2.setItem(i, j, QTableWidgetItem(str(format(float(self.piezoMaterial.epsilon[i, j]),'.1f'))))




    

        # self.lineEdit_v.setText(str(self.materials[item]["velocity"]))
        # self.lineEdit_K2.setText(str(self.materials[item]["K2"]))
        # self.lineEdit_Co.setText(str(self.materials[item]["Co"]))

app = QApplication(sys.argv)
window = MyMainWindow()
window.show()
sys.exit(app.exec())