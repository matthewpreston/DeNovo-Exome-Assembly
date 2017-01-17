#USAGE: python $0 Anolis_Ref_Final_$assembler{,/Summary.csv}

import csv
import glob
import sys

genes = set()
organisms = sorted(glob.glob("%s/*/" % sys.argv[1]))
organisms = [o.split('/')[-2] for o in organisms]
columns = []

for arg in sorted(glob.glob("%s/*/Summary_Results.csv" % sys.argv[1])):
    reader = csv.reader(open(arg, 'r'))
    records = list(reader)[1:-2]
    genes |= set([r[0] for r in records])
    data = {r[0]: float(r[-1]) for r in records}
    columns.append(data)

genes = sorted(list(genes), key=lambda x: x.lower())
orgave = [str(sum(c.values()) / len(genes)) for c in columns]

with open(sys.argv[2], 'w') as o:
    o.write("Gene," + ",".join(organisms) + ",,Gene.Average\n")
    for gene in genes:
        row = []
        for col in columns:
            try:
                row.append(str(col[gene]))
            except KeyError:
                row.append("0")
        geneave = str(sum([float(r) for r in row]) / len(organisms))
        o.write(gene + "," + ",".join(row) + ",," + geneave + "\n")
    o.write("\nAverage," + ",".join(orgave) + "\n")