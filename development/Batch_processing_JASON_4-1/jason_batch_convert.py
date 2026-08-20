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

import os
import datetime
from multiprocessing import Pool, Manager
from multiprocessing.managers import ValueProxy
from pathlib import Path
import argparse
import glob
from colorama import Fore, init

import beautifuljason as bjason

def unique_output_filename(file, extension):
    """Generate a unique output filename based on the input file path and desired extension."""
    unique_base_name = file.replace(os.sep, '_')  # Replace path separators with underscores
    return unique_base_name + '.' + extension

def process_file(execute: bool, file: str, in_dir: str, out_dir: str, jason: bjason.JASON, formats: list[str], rules: str, counter: ValueProxy, total_files: int):
    """Process a single file, converting it to specified formats."""
    out_fnames = [unique_output_filename(file, format) for format in formats]
    if execute:
        global processed_files
        with jason.create_document(os.path.join(in_dir, file), rules=rules) as doc:
            doc.close()
            jason.save(doc, [os.path.join(out_dir, fname) for fname in out_fnames])
    counter.value += 1
    print(f"{counter.value}/{total_files}: {file} => {', '.join(out_fnames)}")

def main():
    init(autoreset=True)  # Initialize colorama

    parser = argparse.ArgumentParser(
        description='Convert files from a specified directory based on extensions or patterns.',
        usage='jason_batch_convert [-h] in_dir out_dir --formats FORMATS [FORMATS...] (--extensions EXTENSIONS [EXTENSIONS ...] | --patterns PATTERNS [PATTERNS ...]) [--rules RULES] [--execute]'
        )
    parser.add_argument('in_dir', help='Root directory containing files to be converted.')
    parser.add_argument('out_dir', help='Directory where the converted files will be saved.')
    parser.add_argument('--formats', nargs='+', required=True, 
                        choices=['jjh5', 'jjj', 'jdx', 'jdf', 'pdf', 'png', 'jpg', 'svg'],
                        help='List of desired output file formats. Multiple formats can be specified.')
    parser.add_argument('--rules', type=str, default="off",
                        help='Path to the rules file or rule library name.')
    parser.add_argument('--execute', action='store_true', 
                        help='Execute the file conversions. If not specified, will perform a dry run, listing the files to be converted. It is recommended to run a dry run first.')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--extensions', nargs='+', help='List of file extensions to convert, e.g., jdf, als, jdx, mol, etc.')
    group.add_argument('--patterns', nargs='+', help='List of file path patterns to convert, e.g., */*.fid/fid, */10/fid, etc. Patterns cannot end with "/*".')

    args = parser.parse_args()

    # Convert in_dir and out_dir to absolute paths
    args.in_dir = os.path.abspath(args.in_dir)
    args.out_dir = os.path.abspath(args.out_dir)

    if not args.execute:
        print(f"{Fore.YELLOW}Dry run mode! No files will be converted.")
    print("Input directory:", args.in_dir)
    print(f"Output directory:", args.out_dir)

    # Initialize Jason
    jason = bjason.JASON()

    # Create output directory
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)

    # Find and process files matching patterns or extensions
    print("Searching for files...")
    files_to_process = []
    if args.patterns:
        for pattern in args.patterns:
            if pattern.endswith("/*"):
                raise ValueError("Invalid pattern. Patterns cannot end with '/*'.")
            pattern_files = glob.glob(os.path.join(args.in_dir, pattern), recursive=True)
            files_to_process.extend([os.path.relpath(f, start=args.in_dir) for f in pattern_files if os.path.isfile(f)])
    else:
        for root, dirs, files in os.walk(args.in_dir):
            for file in files:
                if file.split('.')[-1] in args.extensions:
                    full_path = os.path.join(root, file)
                    if os.path.isfile(full_path):
                        files_to_process.append(os.path.relpath(full_path, start=args.in_dir))

    total_files = len(files_to_process)
    print(f"{total_files} files found in the input directory.")

    manager = Manager()
    counter = manager.Value('i', 0)  # integer counter initialized to 0

    if args.execute:
        start_time = datetime.datetime.now()
        print("Conversion started at:", start_time)
    with Pool() as pool:
        pool.starmap(process_file, [(args.execute, file, args.in_dir, args.out_dir, jason, args.formats, args.rules, counter, total_files) for file in files_to_process])
    if args.execute:
        end_time = datetime.datetime.now()
        print("Conversion finished at:", end_time)
        print("Duration:", end_time - start_time)

if __name__ == '__main__':
    main()

