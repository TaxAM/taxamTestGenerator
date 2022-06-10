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
        '-n', '--pool_name',
        help = 'Pool name',
        type = str,
        action = 'store',
        default = None
    )

    parser.add_argument(
        '-s', '--sample_names',
        help = 'Sample names',
        type = str,
        action = 'store',
        default = None
    )

    parser.add_argument(
        '-nt', '--number_of_taxa_per_level',
        help = 'Number of taxa per level',
        type = str,
        action = 'store',
        default = None
    )

    parser.add_argument(
        '-pt', '--percentage_of_partial_taxa',
        help = 'Percentage of partial taxa',
        type = float,
        action = 'store',
    )

    parser.add_argument(
        '-nr', '--number_of_reads',
        help = 'Number of reads',
        type = int,
        action = 'store',
    )

    parser.add_argument(
        '-nc', '--number_of_contigs',
        help = 'Number of contigs',
        type = int,
        action = 'store',
    )

    parser.add_argument(
        '-pm', '--percentage_of_mapped_reads',
        help = 'Percentage of mapped reads',
        type = float,
        action = 'store',
    )

    parser.add_argument(
        '-tr', '--number_of_taxa_per_read',
        help = 'Number of taxa per read',
        type = int,
        action = 'store',
    )

    parser.add_argument(
        '-tc', '--number_of_taxa_per_contig',
        help = 'Number of taxa per contig',
        type = int,
        action = 'store',
    )

    parser.add_argument(
        '-cr', '--percentage_of_classified_reads',
        help = 'Percentage of classified reads',
        type = float,
        action = 'store',
    )

    parser.add_argument(
        '-cc', '--percentage_of_classified_contigs',
        help = 'Percentage of classified contigs',
        type = float,
        action = 'store',
    )

    parser.add_argument(
        '-mc', '--percentage_of_matched_class',
        help = 'Percentage of matched class',
        type = float,
        action = 'store',
    )

    return parser.parse_args()