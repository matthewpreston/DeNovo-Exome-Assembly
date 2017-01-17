# Adds a prefix string to each FASTA record (can combine multiple FASTA files)
# 
# Written By: Matt Preston (matthew.preston@mail.utoronto.ca)
# Created On: Oct 6, 2016
# Revised On: Never

import argparse
import sys
from Bio import SeqIO

USAGE  = """Usage: python %(prog)s [options] <input.fasta>..."""
DESC   = """Adds a prefix string to each FASTA record (can combine multiple FASTA
            files)"""
EPILOG = """If bugs are found or you wish for more features to be added,
            email: matthew.preston@mail.utoronto.ca"""

def type_prefix(prefix):
    """Tests if prefix is valid (not an empty string)"""
    if prefix == "":
        raise argparse.ArgumentTypeError("Prefix is an empty string!")
    return prefix

def main(): #Gasp! No argument vector anymore!
    #Set up parser with explantory text
    parser = argparse.ArgumentParser(usage=USAGE,
                                     description=DESC,
                                     epilog=EPILOG)
    #Add optional arguments
    parser.add_argument("-o", "--output",
                        default="out.txt",
                        help="Name of the output file (%(default)s)",
                        metavar="FILE")
    parser.add_argument("-p", "--prefix",
                        default="",
                        type=type_prefix,
                        help="Prefix string to add to each record"
                             "(\"%(default)s\")",
                        metavar="STR")
    parser.add_argument("-r", "--replace_str",
                        nargs=2,
                        help="Replaces all first string occurences with the "
                             "second. This step is done before the prefix is "
                             "added.",
                        metavar="STR")
    #Add positional arguments
    parser.add_argument("inputs",
                        nargs='+',
                        help="File containing FASTA records to be prefixed",
                        metavar="FILE")
    #Parse command line (default is sys.argv)
    args = parser.parse_args()
    
    #Parse through records, add prefix, fix description, write to file
    with open(args.output, 'w') as o:
        for input in args.inputs:
            for record in SeqIO.parse(input, "fasta"):
                #Add prefix
                if hasattr(args, "replace_str"):
                    record.id = args.prefix \
                                + record.description.replace(*args.replace_str)
                else:
                    record.id = args.prefix + record.description
                #Fix description b/c SeqIO.write is dumb (see its documentation)
                record.description = ""
                SeqIO.write(record, o, "fasta")

if __name__ == "__main__":
    main()