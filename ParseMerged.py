import sys
from Bio import SeqIO


def main(argv):
    outfile = argv[0]
    records = []
    for arg in argv[1:]:
        records.append(list(SeqIO.parse(arg, "fasta"))[1])
    SeqIO.write(records, outfile, "fasta")


if __name__ == "__main__":
    main(sys.argv[1:])
