# Loading_inputdata.py
from PyQt5.QtWidgets import QFileDialog
import json  # Import the json module
import argparse

def browse_file(ui):
    options = QFileDialog.Options()
    file_path, _ = QFileDialog.getOpenFileName(None, "Select the input file:", "",
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

        if faults:
            # Access the loaded JSON data (example)
            print(f"Loaded JSON data: {len(faults)} faults found.")


# Call your function with the input data
result = process_json_file(input_data)
