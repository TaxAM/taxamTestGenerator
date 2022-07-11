"""That's just a script to convert matix outputs to the right format
"""
# Built-in imports
import glob, os

# Local imports
from utils.local_parsers import local_parsers

args = local_parsers()

files_folder = args.input_folder

files_to_covert = glob.glob(files_folder + '/' + args.file_regex)

if len(files_to_covert):
    for file in files_to_covert:

        output_folder = args.output_folder
        if(not os.path.isdir(output_folder)):
            os.mkdir(output_folder)

        with open(file, 'r') as fr, open(output_folder + '/' + os.path.split(file)[-1], 'w') as fw:
            file_lines = fr.readlines()
            print(file)
            header = file_lines[0]
            matrix = {}
            for line in file_lines[1:]:
                key, value1, value2 = line.replace('\n', '').split('\t')
                key = key.split('; ')[-1].strip()
                values = [int(value1), int(value2)]
                try:
                    matrix[key][0] += values[0]
                    matrix[key][1] +=  values[1]
                except:
                    matrix[key] = values

            matrix_lines = header
            sorted_keys = sorted(matrix.keys())
            for key in sorted_keys:
                matrix_lines += '\t'.join([key, *[str(value) for value in matrix[key]]]) + '\n'
            fw.write(matrix_lines)
            print('-'*80)
else:
    exit('There is no one file to be converted!')