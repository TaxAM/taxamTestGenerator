"""Module to store functions about local parsers arguments."""
import argparse

def local_parsers():
    """Function to create and return all the local parsers arguments.

    Returns
    -------
    argparse
        Object with parsers and its values.
    """
    parser = argparse.ArgumentParser(
        description = 'That s a tool to generate test data for the Taxam project. It creates and counts the taxa.',
        usage = 'python ./taxamTestGenerator <FLAGS>'
    )

    parser.add_argument(
        '-t', '-test',
        help = 'Just a test parser',
        type = str,
        action = 'store',
        default = None
    )