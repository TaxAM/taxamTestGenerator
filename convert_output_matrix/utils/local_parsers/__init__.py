"""Module to store functions about local parsers arguments."""
import argparse

def local_parsers():
    """Function to create and return all the local parsers arguments.

    Returns
    -------
    argparse
        Object with parsers and its values.
    """    
    parser = argparse.ArgumentParser(description='Converts a matrix where each animal is a group with more than one animal for a matrix where each animal will be just a taxon level.', usage='',)

    parser.add_argument(
        '-i', '--input_folder',
        help = 'Folder where are the files to be converted.',
        type = str,
        action = 'store',
        default = None,
    )

    parser.add_argument(
        '-f', '--file_regex',
        help = 'Regex used in the folder to find the files to be converted.',
        type = str,
        action = 'store',
        default = '*.taxam',
    )

    parser.add_argument(
        '-o', '--output_folder',
        help = 'Folder where store converted files.',
        type = str,
        action = 'store',
        default = r'./converted_files/',
    )

    return parser.parse_args()



