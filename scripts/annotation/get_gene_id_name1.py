import re
import os
import sys

infile, outfile = sys.argv[1:]

f1 = open(infile)
id_name_dct = {}
for line in f1:
    if line.startswith("#"):
        continue
    line = line.split("\t")[8]
    m = re.findall(r"gene_id\s+\"([^;]+)\";", line)
    n = re.findall(r"gene_id\s+\"([^;]+)\";", line)
    if m != [] and n !=[]:
        gene_id = m[0]
        gene_name = n[0]
        if gene_id not in id_name_dct:
            id_name_dct[gene_id] = gene_name
    if m != [] and n ==[]:
        gene_id =  m[0]
        gene_name = m[0]
        if gene_id not in id_name_dct:
            id_name_dct[gene_id] = gene_name
f1.close()

f2 = open(outfile, "w")
f2.write("gene_id" + "\t"+ "gene_name" + "\n")
for k, v in id_name_dct.items():
    line = k + "\t" + v + "\n"
    f2.write(line)
f2.close()
