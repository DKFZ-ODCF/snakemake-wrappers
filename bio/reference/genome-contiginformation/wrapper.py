__author__ = "Philip R. Kensche, Fritjof Lammers"
__copyright__ = "Copyright 2020, DKFZ"
__email__ = "f.lammers@dkfz.de"
__license__ = "MIT"


# Read in a FASTA file and return a TSV with one entry per FASTA entry and three columns (entryId, totalLength, lengthOnlyATCG).
# This is the information needed for the "stats" file.

import sys
import os.path
from Bio import SeqIO
from collections import Counter

fasta_in = snakemake.input[0]
tab_out = snakemake.output[0]

if not fasta_in:
    sys.exit("Please supply FASTA file, and only FASTA file! Produce a headerless TSV with columns (entryId, totalLength, lenthOnlyATCG.")
if fasta_in:
    if not os.path.isfile(fasta_in):
        sys.exit("Supplied FASTA file cannot be found!")


in_fh = open(fasta_in, "r")
out_fh = open(tab_out, "w")

print("Reading sequence data from " + fasta_in, file=sys.stderr)

for record in SeqIO.parse(in_fh, "fasta"):
    print("Processing entry '{}' ...".format(record.id), file=sys.stderr)
    counter = Counter(record.seq)
    allBases = sum(counter.values())
    onlyATCG = sum(list(map(lambda c: counter[c], ["A", "a", "T", "t", "C", "c", "G", "g"])))
    print("\t".join([record.id, str(allBases), str(onlyATCG)]), file=out_fh)

in_fh.close()
out_fh.close()
