# FaultQuake_functions.py

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on September 2023

@authors: Nasrin Tavakolizadeh, Hamzeh Mohammadigheymasi
"""

# This module is the main module of the FaultQuake tool to calculate the seismic activity rates of faults.




from src.SeismicActivityRate import momentbudget, sactivityrate # Import the function
import sys, argparse, json
global faults
from PyQt5.QtWidgets import QFileDialog
from PyQt5 import QtCore, QtGui, QtWidgets

class InfoWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Information")
        self.setGeometry(100, 100, 600, 400) # Adjusted size for better visibility
        layout = QtWidgets.QVBoxLayout()
        info_text = """
        {
            "MFF1": {
                "ScR": "WC94-R",
                "year_for_calculations": 2023,
                "Length": 37,
                "Dip": 40,
                "Seismogenic_Thickness": 20,
                "SRmin": 2.88,
                "SRmax": 4.32,
                "Mobs": 5.7,
                "sdMobs": 0.05,
                "Last_eq_time": 1987,
                "SCC": 20.5,
                "ShearModulus": "NaN",
                "StrainDrop": 3,
                "Mmin": 5.5,
                "b-value": 0.9
            },
        }
        """
        # Use a scrollable area if the text is too large
        scroll_area = QtWidgets.QScrollArea()
        self.label = QtWidgets.QLabel(info_text)
        self.label.setWordWrap(True)
        scroll_area.setWidget(self.label)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        self.setLayout(layout)


class ScaleWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Information")
        self.setGeometry(100, 100, 600, 400)  # Adjusted size for better visibility
        layout = QtWidgets.QVBoxLayout()

        # Replace the existing info_text with the new content
        info_text = """
        set ScR as follows

        Wells and Coppersmith (1984) relationships:
        WC94-N - normal faults
        WC94-R - reverse faults
        WC94-S - strike slip faults
        WC94-A - all the kinematics

        Leonard (2010) relationships:
        Le10-D - dip slip faults
        Le10-S - strike slip faults
        Le10-SCR - stable continental regions

        Volcanic context relationships (Azzaro et al., 2015; Villamor et al.,
        2001):
        Volc - all the kinematics
        """
        # Use a scrollable area if the text is too large
        scroll_area = QtWidgets.QScrollArea()
        self.label = QtWidgets.QLabel(info_text)
        self.label.setWordWrap(True)
        scroll_area.setWidget(self.label)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        self.setLayout(layout)
class OutWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Information")
        self.setGeometry(100, 100, 600, 400) # Adjusted size for better visibility
        layout = QtWidgets.QVBoxLayout()
        info_text = """
<?xml version="1.0" encoding="utf-8"?>
<nrml xmlns="http://openquake.org/xmlns/nrml/0.4" xmlns:gml="http://www.opengis.net/gml">
    <sourceModel name="FFF1">
        <simpleFaultSource id="1" name="Simple Fault Source" tectonicRegion="Active Shallow Crust">
            <simpleFaultGeometry>
                <gml:LineString>
                    <gml:posList>
                        <!-- Fill in coordinates here -->
                    </gml:posList>
                </gml:LineString>
                <dip>40</dip>
                <upperSeismoDepth>-1.600000e+01</upperSeismoDepth>
                <lowerSeismoDepth>2.000000e+01</lowerSeismoDepth>
            </simpleFaultGeometry>
            <magScaleRel>WC94-R</magScaleRel>
            <ruptAspectRatio>2.0000000E+00</ruptAspectRatio>
                <occurRates>2.325063e-01 1.889882e-01 1.536154e-01 1.248633e-01 1.014927e-01 8.249635e-02 6.705555e-02 5.450480e-02 4.430316e-02 3.601096e-02 2.927081e-02 2.379221e-02 1.933903e-02</occurRates>
            </incrementalMFD>
            <rake>9.0000000E+01</rake>
        </simpleFaultSource>
    </sourceModel>
</nrml>
        """
        # Use a scrollable area if the text is too large
        scroll_area = QtWidgets.QScrollArea()
        self.label = QtWidgets.QLabel(info_text)
        self.label.setWordWrap(True)
        scroll_area.setWidget(self.label)
        scroll_area.setWidgetResizable(True)
        layout.addWidget(scroll_area)
        self.setLayout(layout)


def browse_file(ui, self=None):
    options = QFileDialog.Options()
    file_path, _ = QFileDialog.getOpenFileName(None, "Select the input file:", "Choose one:",
                                               "JSON files (*.json);;All files (*)", options=options)
    if file_path:
        print("Selected file:", file_path)  # Print the selected file path
        # Store the selected file path in the class variable
        ui.input_file_var = file_path

        # Load the selected JSON file
        try:
            with open(file_path, 'r') as f:
                faults = json.load(f)
                ui.faults = faults

            if faults:
                # Access the loaded JSON data (example)
                print(f"Loaded JSON data: {len(faults)} faults found.")
        except Exception as e:
            print("Error loading JSON file:", e)


class Ui_Frame(object):
    def setupUi(self, Frame):
        Frame.setObjectName("FaultQuake")
        Frame.resize(903, 758)
        self.horizontalLayout = QtWidgets.QHBoxLayout(Frame)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.frame = QtWidgets.QFrame(Frame)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(233, 185, 110))
        brush.setStyle(QtCore.Qt.SolidPattern)
        self.frame.setPalette(palette)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.pushButton = QtWidgets.QPushButton(self.frame)
        self.pushButton.clicked.connect(lambda: browse_file(ui))
        self.pushButton.setGeometry(QtCore.QRect(140, 80, 141, 41))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        self.pushButton.setPalette(palette)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(17)
        font.setBold(False)
        font.setWeight(50)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtWidgets.QPushButton(self.frame)
        self.pushButton_2.setGeometry(QtCore.QRect(120, 500, 201, 71))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(204, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Button, brush)
        self.pushButton_2.setPalette(palette)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(24)
        font.setBold(True)
        font.setWeight(75)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        # Use a style sheet to set the background color of the Run button
        self.pushButton_2.setStyleSheet("background-color: navy; color: white;")
        # Connect the "RUN" button to call_my_function
        self.pushButton_2.clicked.connect(lambda: self.SeismicActivityRate(self.faults, self.mfdo))
        # MFD Options combe box
        self.comboBox = QtWidgets.QComboBox(self.frame)
        self.comboBox.setGeometry(QtCore.QRect(30, 250, 391, 41))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(188, 188, 202))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Button, brush)
        self.comboBox.setPalette(palette)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(16)
        self.comboBox.setFont(font)
        self.comboBox.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.comboBox.setObjectName("comboBox")
        self.comboBox.currentIndexChanged.connect(self.set_mfdo)  # Connect the combobox to the set_mfdo method
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.label = QtWidgets.QLabel(self.frame)
        self.label.setGeometry(QtCore.QRect(60, 190, 291, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(16)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.label.setTextFormat(QtCore.Qt.AutoText)
        self.label.setIndent(-3)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.frame)
        self.label_2.setGeometry(QtCore.QRect(120, 30, 181, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(16)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_2.setFont(font)
        self.label_2.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.label_2.setTextFormat(QtCore.Qt.AutoText)
        self.label_2.setIndent(-3)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.frame)
        self.label_3.setGeometry(QtCore.QRect(100, 350, 241, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(16)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.label_3.setTextFormat(QtCore.Qt.AutoText)
        self.label_3.setIndent(-3)
        self.label_3.setObjectName("label_3")
        self.textEdit = QtWidgets.QTextEdit(self.frame)
        self.textEdit.setGeometry(QtCore.QRect(130, 400, 181, 41))
        self.textEdit.setPalette(palette)
        self.textEdit.setObjectName("textEdit")
        self.horizontalLayout.addWidget(self.frame)
        self.frame_2 = QtWidgets.QFrame(Frame)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(211, 215, 207))
        brush.setStyle(QtCore.Qt.SolidPattern)
        self.frame_2.setPalette(palette)
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.label_4 = QtWidgets.QLabel(self.frame_2)
        self.label_4.setGeometry(QtCore.QRect(10, 430, 81, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(16)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.label_4.setFont(font)
        self.label_4.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.label_4.setTextFormat(QtCore.Qt.AutoText)
        self.label_4.setIndent(-3)
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.frame_2)
        self.label_5.setGeometry(QtCore.QRect(10, 350, 261, 21))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        font.setBold(True)
        font.setItalic(False)
        font.setWeight(75)
        self.label_5.setFont(font)
        self.label_5.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.label_5.setTextFormat(QtCore.Qt.AutoText)
        self.label_5.setIndent(-3)
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.frame_2)
        self.label_6.setGeometry(QtCore.QRect(10, 550, 91, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(13)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_6.setFont(font)
        self.label_6.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.label_6.setTextFormat(QtCore.Qt.AutoText)
        self.label_6.setIndent(-3)
        self.label_6.setObjectName("label_6")
        self.textEdit_3 = QtWidgets.QTextEdit(self.frame_2)
        self.textEdit_3.setGeometry(QtCore.QRect(110, 390, 201, 31))
        self.textEdit_3.setPalette(palette)
        self.textEdit_3.setObjectName("textEdit_3")

        self.textEdit_4 = QtWidgets.QTextEdit(self.frame_2)
        self.textEdit_4.setGeometry(QtCore.QRect(110, 470, 201, 31))
        self.textEdit_4.setPalette(palette)
        self.textEdit_4.setObjectName("textEdit_4")


        # Descriptions pannel
        self.pushButton_3 = QtWidgets.QPushButton(self.frame_2)
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.clicked.connect(self.open_info_window)
        self.pushButton_3.setGeometry(QtCore.QRect(100, 590, 231, 41))
        font = QtGui.QFont("Times New Roman", 14)
        font.setFamily("Times New Roman")
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_3.setFont(font)
        self.pushButton_3.setObjectName("pushButton_3")

        self.pushButton_5 = QtWidgets.QPushButton(self.frame_2)
        self.pushButton_5.setGeometry(QtCore.QRect(100, 640, 231, 41))
        font = QtGui.QFont("Times New Roman", 14)
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_5.setFont(font)
        self.pushButton_5.setObjectName("pushButton_5")
        self.pushButton_5.clicked.connect(self.open_scale_window)


        self.pushButton_6 = QtWidgets.QPushButton(self.frame_2)
        self.pushButton_6.setGeometry(QtCore.QRect(100, 690, 231, 41))
        font = QtGui.QFont("Times New Roman", 14)
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_6.setFont(font)
        self.pushButton_6.setObjectName("pushButton_6")
        self.pushButton_6.clicked.connect(self.open_out_info_window)


        self.frame_3 = QtWidgets.QFrame(self.frame_2)
        self.frame_3.setGeometry(QtCore.QRect(30, 60, 371, 171))
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")


        self.label_8 = QtWidgets.QLabel(self.frame_3)
        self.label_8.setGeometry(QtCore.QRect(10, 10, 171, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(15)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")


        # Zeta

        self.textEdit_5 = QtWidgets.QTextEdit(self.frame_3)
        self.textEdit_5.setGeometry(QtCore.QRect(160, 110, 101, 31))
        self.textEdit_5.setPalette(palette)
        self.textEdit_5.setObjectName("textEdit_5")



        # Khi:

        self.textEdit_6 = QtWidgets.QTextEdit(self.frame_3)
        self.textEdit_6.setGeometry(QtCore.QRect(160, 60, 101, 31))
        self.textEdit_6.setPalette(palette)
        self.textEdit_6.setObjectName("textEdit_6")



        # Siggma

        self.textEdit_7 = QtWidgets.QTextEdit(self.frame_2)
        self.textEdit_7.setGeometry(QtCore.QRect(110, 290, 201, 31))
        self.textEdit_7.setPalette(palette)
        self.textEdit_7.setObjectName("textEdit_7")



        self.label_9 = QtWidgets.QLabel(self.frame_3)
        self.label_9.setGeometry(QtCore.QRect(90, 60, 51, 31))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")


        self.label_10 = QtWidgets.QLabel(self.frame_3)
        self.label_10.setGeometry(QtCore.QRect(90, 110, 41, 21))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")


        self.label_7 = QtWidgets.QLabel(self.frame_2)
        self.label_7.setGeometry(QtCore.QRect(130, 10, 161, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(17)
        font.setBold(True)
        font.setWeight(75)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")


        self.label_11 = QtWidgets.QLabel(self.frame_2)
        self.label_11.setGeometry(QtCore.QRect(10, 260, 251, 21))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")




        self.horizontalLayout.addWidget(self.frame_2)
        self.retranslateUi(Frame)
        QtCore.QMetaObject.connectSlotsByName(Frame)



    def SeismicActivityRate(self, faults, mfdo):
        ProjFol = self.textEdit.toPlainText()
        PTI = self.textEdit_3.toPlainText()
        MBS = self.textEdit_4.toPlainText()
        Khi = self.textEdit_5.toPlainText()
        Zeta = self.textEdit_6.toPlainText()
        Siggma = self.textEdit_7.toPlainText()

        # Set default values if empty
        if not PTI:
            PTI = "50"
            print('Window of observation: Where necessary you are using default value 50 years')
        if not MBS:
            MBS = "0.1"
            print('binstep: Where necessary you are using default value 0.1')
        try:
            PTI = float(PTI)
            MBS = float(MBS)
        except ValueError:
            raise ValueError(
                "Consider inputting a proper number for 'Probability Time Interval', and 'Magnitude Bin Size'")

        if not Zeta:
            Zeta = "0.5"
            print('Magnitude difference: Where necessary you are using default value 0.5')
        if not Khi:
            Khi = "0.2"
            print('Multiplication Factor: Where necessary you are using default value 0.2')

        if not Siggma:
            Siggma = "0.3"
            print('Standard Deviation: Where necessary you are using default value 0.3')
        try:
            Khi = float(Khi)
            Zeta = float(Zeta)
            Siggma = float(Siggma)

        except ValueError:
            raise ValueError(
                "Consider inputting a proper number for 'Probability Time Interval', and 'Magnitude Bin Size'")




        faults_u = momentbudget(faults, Zeta, Khi, Siggma, ProjFol='output_files', logical_nan='NAN, "",NaN',  logical_nan_sdmag='NAN, "",NaN')
        if ProjFol == '':
            ProjFol == 'output_files'
        sactivityrate(faults_u, mfdo, PTI, MBS, ProjFol='output_files')

    def set_mfdo(self, index):
        options = [
            "Magnitude Frequency Distribution options",
            "Truncated Gutenberg Richter",
            "Characteristic Gaussian",

        ]
        mfdo = options[index]  # Get the selected option
        self.mfdo = mfdo  # Store the selected value for later use in the SeismicActivityRate method

    def open_info_window(self):
        self.info_window = InfoWindow()
        self.info_window.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint)  # Set the window to stay on top
        self.info_window.show()

    def open_out_info_window(self):
        self.info_window = OutWindow()
        self.info_window.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint)  # Set the window to stay on top
        self.info_window.show()

    def open_scale_window(self):
        self.scale_window = ScaleWindow()  # You can create a separate InfoWindow for the scale window if needed
        self.scale_window.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint)  # Set the window to stay on top
        self.scale_window.show()



    def retranslateUi(self, Frame):
        _translate = QtCore.QCoreApplication.translate
        Frame.setWindowTitle(_translate("FaultQuake", "FaultQuake"))
        self.pushButton.setText(_translate("FaultQuake", "Browse"))
        self.pushButton_2.setText(_translate("FaultQuake", "RUN"))
        self.comboBox.setItemText(0, _translate("FaultQuake", "Magnitude Frequency Distribution options"))
        self.comboBox.setItemText(1, _translate("FaultQuake", "Truncated Gutenberg Richter"))
        self.comboBox.setItemText(2, _translate("FaultQuake", "Characteristic Gaussian"))
        self.label.setText(_translate("FaultQuake", "Activity Rate calculation options:"))
        self.label_2.setText(_translate("FaultQuake", "Select the input file:"))
        self.label_3.setText(_translate("FaultQuake", "Write the output file name:"))
        self.label_4.setText(_translate("FaultQuake", "bin size:"))
        self.label_5.setText(_translate("FaultQuake", "Probability Window (years):"))
        self.label_6.setText(_translate("FaultQuake", "Discriptions:"))
        self.pushButton_3.setText(_translate("FaultQuake", "FaultQuake Input file format"))
        self.pushButton_5.setText(_translate("FaultQuake", "Scale-Relationship selection"))
        self.pushButton_6.setText(_translate("FaultQuake", "FaultQuke output file format"))
        self.label_8.setText(_translate("FaultQuake", "OVCW Parameters:"))
        self.label_9.setText(_translate("FaultQuake", "<html>&#950; =</html>"))  # Greek letter xi: ξ
        self.label_10.setText(_translate("FaultQuake", "<html>&#958; =</html>"))  # Greek letter zeta: ζ
        self.label_7.setText(_translate("FaultQuake", "Optional Inputs:"))
        self.label_11.setText(_translate("FaultQuake", "M<sub>w</sub>(M<sub>0</sub>) Standard Deviation:"))
        self.textEdit_3.setPlaceholderText(_translate("FaultQuake", "50"))
        self.textEdit_4.setPlaceholderText(_translate("FaultQuake", "0.1"))
        self.textEdit_5.setPlaceholderText(_translate("FaultQuake", "0.2"))
        self.textEdit_6.setPlaceholderText(_translate("FaultQuake", "0.5"))
        self.textEdit_7.setPlaceholderText(_translate("FaultQuake", "0.3"))


if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    Frame = QtWidgets.QFrame()
    ui = Ui_Frame()
    ui.setupUi(Frame)
    Frame.show()
    Frame.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint)
    Frame.show()
    sys.exit(app.exec_())
