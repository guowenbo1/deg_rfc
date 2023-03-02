import sys
import os
import re

infile = sys.argv[1]
outfile = sys.argv[2]

ext = os.path.splitext(infile)[-1]
gene_id_pattern = ""
print(ext)
if ext in ".gff":
	gene_id_pattern = "ID=(.*?);"
elif ext == ".gtf":
	gene_id_pattern = "gene_id \"(.*?)\";"


with open(infile, "r") as in_f_1:
	anno = [i for i in in_f_1.readlines() if not i.startswith("#") and i.strip().split("\t")[2] == "gene"]


with open(outfile, "w") as out_f_1:
	out_f_1.write("gene_id\tchr\tStart\tEnd\tStrand\n")

	for i in anno:
		i_lst = i.strip().split("\t")
		gene_id = re.findall(gene_id_pattern, i_lst[8])[0]
		chrosome = i_lst[0].replace("Chr", "chr")
		if "chr" not in chrosome:
			chrosome = "chr" + chrosome
		Start = i_lst[3]
		End = i_lst[4]
		Strand = i_lst[6]

		out_f_1.writelines("\t".join(list(map(str, [gene_id, chrosome, Start, End, Strand]))) + "\n")
