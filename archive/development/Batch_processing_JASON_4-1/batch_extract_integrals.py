##
##-----------------------------------------------------------------------------
##
## Copyright (c) 2024 JEOL Ltd.
## 1-2 Musashino 3-Chome
## Akishima Tokyo 196-8558 Japan
##
## This software is provided under the MIT License. For full license information,
## see the LICENSE file in the project root or visit https://opensource.org/licenses/MIT
##
##++---------------------------------------------------------------------------
##
## ModuleName : BeautifulJASON
## ModuleType : Python API for JASON desktop application and JJH5 documents
## Purpose : Automate processing, analysis, and report generation with JASON
## Author : Nikolay Larin
## Language : Python
##
####---------------------------------------------------------------------------
##

import argparse
import os
import csv
import io
import beautifuljason as bjason

def escape_newlines_for_csv(value):
    """Escapes newline characters in a string for CSV compatibility."""
    value = str(bjason.utils.ensure_str(value))
    # Replace newline characters with literal \n to ensure CSV compatibility
    return value.replace('\n', '\\n').replace('\r', '\\r')
    
def extract_integrals(jjh5_file_path: str, csv_writer: csv.writer, params: list[str]):
    """Extracts integrals from a JJH5 file and writes them to a CSV file using csv.writer."""
    try:
        base_name = os.path.basename(jjh5_file_path)
        file_name_without_extension, _ = os.path.splitext(base_name)
        
        with bjason.Document(jjh5_file_path, mode="r") as doc:
            spec = doc.nmr_data[0]
            row = [file_name_without_extension]  # Start row with file name
            
            for param in params:
                param_parts = param.split('/')
                param_group = param_parts[0]
                param_name = param_parts[1]
                param_index = None
                if '[' in param_name:
                    param_name, param_index = param_name.split('[')
                    param_index = int(param_index[:-1])

                if param_group == "jason_parameters":
                    param_value = spec.spec_info.get_param(param_name)
                else:
                    param_value = spec.raw_data.spec_info.get_orig_param(param_group, param_name)

                if param_index is not None and param_value is not None:
                    if hasattr(param_value, '__getitem__'):
                        param_value = param_value[param_index]

                row.append(escape_newlines_for_csv(param_value))
                
            row.extend(integral.value_hz for integral in spec.multiplets)
            csv_writer.writerow(row)
    except Exception as e:
        print(f"Error processing file '{jjh5_file_path}': {e}")

def iterate_jjh5_files(jjh5_path, csv_path, params):
    if not os.path.isdir(jjh5_path):
        print(f"The directory '{jjh5_path}' does not exist.")
        return
    if not os.listdir(jjh5_path):
        print("The directory is empty.")
        return

    try:
        with open(csv_path, "w", newline='', encoding="utf-8") as csv_file:
            csv_writer = csv.writer(csv_file, quoting=csv.QUOTE_MINIMAL)
            found_files = False
            for file_name in os.listdir(jjh5_path):
                if file_name.endswith(".jjh5"):
                    found_files = True
                    print(f"Processing file: {file_name}")
                    extract_integrals(os.path.join(jjh5_path, file_name), csv_writer, params)
            if not found_files:
                print("No .jjh5 files found in the directory.")
    except IOError as e:
        print(f"Failed to write to CSV file {csv_path}: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Iterates through .jjh5 files in a directory and extracts integrals to a CSV file.",
        usage="batch_extract_integrals [-h] jjh5_path csv_path [-p PARAMETERS ...]",
        epilog="example: batch_extract_integrals path/to/jjh5_files path/to/output.csv -p parameters/ACTUAL_START_TIME jason_parameters/SpectrometerFrequencies[0]"
        )
    parser.add_argument("jjh5_path", type=str, help="The path to the directory containing .jjh5 files")
    parser.add_argument("csv_path", type=str, help="The path to the output CSV file. If the file exists, it will be overwritten.")
    parser.add_argument("-p", "--parameters", type=str, nargs='+', help="List parameters to extract, each formatted as 'group/parameter' or 'group/parameter[index]'. Ensure to match the exact case used in the JASON GUI. Example: parameters/ACTUAL_START_TIME jason_parameters/SpectrometerFrequencies[0]")
    
    args = parser.parse_args()
    output_dir = os.path.dirname(args.csv_path)
    if not output_dir:
        output_dir = '.'  # If no directory is provided, use the current directory
    if not os.path.exists(output_dir):
        print(f"Output directory '{output_dir}' does not exist.")
        return
    
    iterate_jjh5_files(args.jjh5_path, args.csv_path, args.parameters if args.parameters else [])

if __name__ == "__main__":
    main()
