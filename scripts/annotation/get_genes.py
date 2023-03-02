import sys
import os
import re

infile, outfile = sys.argv[1:]

with open(infile, 'r') as fa:
    sequences = {}
    for line in fa:
        seqs = []
        ids = []
        if line.startswith(">"):
            name = line.rstrip("\n")
            # gene = re.match(r">\S+\s+gene=(\S+).*", name).groups()[0]
            m = re.match(r">\S+\s+gene=([\S\s]+)(CDS).*", name)
            if m is not None:
                gene = m.groups()[0]
            else:
                m = re.match(r">\S+\s+gene=(.*)", name)
                if m is not None:
                     gene = m.groups()[0]
                else:
                    continue
            seq = ""
        else:
            seq = seq + line.rstrip("\n")
        seqs.append(seq)
        if gene not in sequences:
            sequences[gene] = seqs
        else:
            sequences[gene] += seqs

maxseq = {}
for k, v in sequences.items():
    seq = max(v, key=len)
    maxseq[k] = seq

with open(outfile, "w") as f_obj:
    for k, v in maxseq.items():
        f_obj.write(">" + k + "\n")
        f_obj.write(v+"\n")
