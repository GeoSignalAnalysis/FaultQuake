from src.SeismicActivityRate import momentbudget, sactivityrate # Import the function
import sys
global faults
from PyQt5.QtWidgets import QFileDialog
from PyQt5 import QtCore, QtGui, QtWidgets
import argparse
import json


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

def browse_file(ui, self=None):
    options = QFileDialog.Options()
    file_path, _ = QFileDialog.getOpenFileName(None, "Select the input file:", "Choose one:",
                                               "JSON files (*.json);;All files (*)", options=options)
    if file_path:
        # Store the selected file path in the class variable
        ui.input_file_var = file_path
        parser = argparse.ArgumentParser(description='read_json_file')
        parser.add_argument('--config-file', dest='config_file', type=str, help='Configuration file path',
                            default='./Input_data/fault_parameters.json')
        args = parser.parse_args()
        with open(args.config_file, 'r') as f:
            faults = json.load(f)
            ui.faults = faults
        if faults:
            # Access the loaded JSON data (example)
            print(f"Loaded JSON data: {len(faults)} faults found.")

class Ui_Frame(object):
    def setupUi(self, Frame):
        Frame.setObjectName("Frame")
        Frame.resize(869, 662)
        self.horizontalLayout = QtWidgets.QHBoxLayout(Frame)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.frame = QtWidgets.QFrame(Frame)
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(200, 185, 110))
        brush.setStyle(QtCore.Qt.SolidPattern)
        self.frame.setPalette(palette)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.pushButton = QtWidgets.QPushButton(self.frame)
        self.pushButton.clicked.connect(lambda: browse_file(ui))
        self.pushButton.setGeometry(QtCore.QRect(140, 50, 141, 41))
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
        # Define class variables to store the loaded JSON data
        self.json_data = None
        self.input_file_var = ""
        # Button 2 (RUN)
        self.pushButton_2 = QtWidgets.QPushButton(self.frame)
        self.pushButton_2.setGeometry(QtCore.QRect(100, 380, 231, 41))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 220))
        brush.setStyle(QtCore.Qt.SolidPattern)
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(22)
        font.setBold(True)
        font.setWeight(75)
        self.pushButton_2.setFont(font)
        self.pushButton_2.setObjectName("pushButton_2")
        # Use a style sheet to set the background color of the Run button
        self.pushButton_2.setStyleSheet("background-color: #0000DC; color: black;")
        # Connect the "RUN" button to call_my_function
        self.pushButton_2.clicked.connect(lambda: self.SeismicActivityRate(self.faults, self.mfdo))
        # MFD Options combe box
        self.comboBox = QtWidgets.QComboBox(self.frame)
        self.comboBox.setGeometry(QtCore.QRect(20, 170, 391, 41))
        palette = QtGui.QPalette()
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.WindowText, brush)
        brush = QtGui.QBrush(QtGui.QColor(188, 188, 202))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.Button, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.BrightText, brush)
        brush = QtGui.QBrush(QtGui.QColor(0, 0, 0))
        brush.setStyle(QtCore.Qt.SolidPattern)
        palette.setBrush(QtGui.QPalette.Active, QtGui.QPalette.ButtonText, brush)
        brush = QtGui.QBrush(QtGui.QColor(255, 255, 255))
        brush.setStyle(QtCore.Qt.SolidPattern)
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
        self.label.setGeometry(QtCore.QRect(70, 120, 321, 41))
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
        self.label_2.setGeometry(QtCore.QRect(120, 0, 181, 41))
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
        self.label_3.setGeometry(QtCore.QRect(100, 260, 241, 41))
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
        self.textEdit.setGeometry(QtCore.QRect(120, 310, 181, 41))
        self.textEdit.setPalette(palette)
        self.textEdit.setObjectName("textEdit")
        self.horizontalLayout.addWidget(self.frame)
        self.frame_2 = QtWidgets.QFrame(Frame)
        self.frame_2.setPalette(palette)
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.label_4 = QtWidgets.QLabel(self.frame_2)
        self.label_4.setGeometry(QtCore.QRect(10, 110, 151, 41))
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
        self.label_5.setGeometry(QtCore.QRect(10, 10, 311, 21))
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
        self.label_6.setGeometry(QtCore.QRect(10, 420, 121, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(16)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_6.setFont(font)
        self.label_6.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
        self.label_6.setTextFormat(QtCore.Qt.AutoText)
        self.label_6.setIndent(-3)
        self.label_6.setObjectName("label_6")
        self.textEdit_3 = QtWidgets.QTextEdit(self.frame_2)
        self.textEdit_3.setGeometry(QtCore.QRect(10, 40, 351, 31))
        self.textEdit_3.setPalette(palette)
        self.textEdit_3.setObjectName("textEdit_3")
        self.textEdit_4 = QtWidgets.QTextEdit(self.frame_2)
        self.textEdit_4.setGeometry(QtCore.QRect(10, 150, 351, 31))
        self.textEdit_4.setPalette(palette)
        self.textEdit_4.setObjectName("textEdit_4")
        self.pushButton_3 = QtWidgets.QPushButton(self.frame_2)
        self.pushButton_3.setGeometry(QtCore.QRect(60, 560, 291, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(17)
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_3.setFont(font)
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_5 = QtWidgets.QPushButton(self.frame_2)
        self.pushButton_5.setGeometry(QtCore.QRect(60, 510, 291, 41))
        font = QtGui.QFont("Times New Roman", 17)
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_5.setFont(font)
        self.pushButton_5.setObjectName("pushButton_5")
        self.pushButton_6 = QtWidgets.QPushButton(self.frame_2)
        self.pushButton_6.setGeometry(QtCore.QRect(60, 460, 291, 41))
        font = QtGui.QFont()
        font.setFamily("Times New Roman")
        font.setPointSize(17)
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_6.setFont(font)
        self.pushButton_6.setObjectName("pushButton_6")
        self.pushButton_6.clicked.connect(self.open_info_window)
        self.horizontalLayout.addWidget(self.frame_2)
        self.retranslateUi(Frame)
        QtCore.QMetaObject.connectSlotsByName(Frame)

    def SeismicActivityRate(self, faults, mfdo):
        ProjFol = self.textEdit.toPlainText()
        PTI = self.textEdit_3.toPlainText()
        MBS = self.textEdit_4.toPlainText()

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
        faults_u = momentbudget(faults,ProjFol='output_files')
        if ProjFol == '':
            ProjFol == 'output_files'

        sactivityrate(faults_u, mfdo, PTI, MBS, ProjFol='output_files')
    def set_mfdo(self, index):
        options = [
            "Magnitude Frequency Distribution options",
            "Characteristic Gaussian",
            "Truncated Gutenberg Richter",
        ]
        mfdo = options[index]  # Get the selected option
        self.mfdo = mfdo  # Store the selected value for later use in the SeismicActivityRate method
    def open_info_window(self):
        self.info_window = InfoWindow()
        self.info_window.show()
    def retranslateUi(self, Frame):
        _translate = QtCore.QCoreApplication.translate
        Frame.setWindowTitle(_translate("FaultQuake", "FaultQuake"))
        self.pushButton.setText(_translate("FaultQuake", "Browse"))
        self.pushButton_2.setText(_translate("FaultQuake", "RUN"))
        self.comboBox.setItemText(0, _translate("FaultQuake", "Choose one:"))
        self.comboBox.setItemText(1, _translate("FaultQuake", "Characteristic Gaussian"))
        self.comboBox.setItemText(2, _translate("FaultQuake", "Truncated Gutenberg Richter"))
        self.label.setText(_translate("FaultQuake", "Magnitude Frequency Distribution:"))
        self.label_2.setText(_translate("FaultQuake", "Select the input file:"))
        self.label_3.setText(_translate("FaultQuake", "Project Folder name:"))
        self.label_4.setText(_translate("FaultQuake", "Magnitude Bin:"))
        self.label_5.setText(_translate("FaultQuake", "Probability Window Interval (years):"))
        self.label_6.setText(_translate("FaultQuake", "Discriptions:"))
        self.pushButton_3.setText(_translate("FaultQuake", "FaultQuake Output file format"))
        self.pushButton_5.setText(_translate("FaultQuake", "Scale-Relationship selection"))
        self.pushButton_6.setText(_translate("FaultQuake", "FaultQuke Input file format"))
if __name__ == "__main__":
    app = QtWidgets.QApplication([])
    Frame = QtWidgets.QFrame()
    ui = Ui_Frame()
    ui.setupUi(Frame)
    Frame.show()
    Frame.setWindowFlag(QtCore.Qt.WindowStaysOnTopHint)
    Frame.show()
    sys.exit(app.exec_())