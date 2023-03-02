import re
import os
import sys

infile, outfile = sys.argv[1:]

name = os.path.basename(infile)
f1 = open(infile)
id_name_dct = {}
h1 = re.search(".gtf", name)
h2 = re.search(".gff3", name)

def get_value(a,b):
    for line in f1:
        if line.startswith("#"):
            continue
        feat = line.split("\t")[2]
        line = line.split("\t")[8]
        h4 = re.search("gene", feat)
        if h4:    
            m = re.findall(r'%s'%a, line)
        # n = re.findall(r'%s'%b, line)
            if m != []:
                n = re.findall(r'%s'%b, line)
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
            if m == []:
                m = re.findall(r"ID=gene-([^;]+);", line)
                n = re.findall(r"Name=([^;]+);", line)
                if m != [] and n !=[]:
                    gene_id = m[0]
                    gene_name = n[0]
                    if gene_id not in id_name_dct:
                        id_name_dct[gene_id] = gene_name
    f1.close()
    return (id_name_dct)

if h1:
    a = 'gene_id\s+\"([^;]+)\";'
    b = 'gene_name\s+\"([^;]+)\";'
    id_name_dct = get_value(a,b)
# if h1:
#     for line in f1:
#         if line.startswith("#"):
#             continue
#         line = line.split("\t")[8]
#         m = re.findall(r"gene_id\s+\"([^;]+)\";", line)
#         n = re.findall(r"gene_name\s+\"([^;]+)\";", line)
#         if m != [] and n !=[]:
#             gene_id = m[0]
#             gene_name = n[0]
#             if gene_id not in id_name_dct:
#                 id_name_dct[gene_id] = gene_name
#         if m != [] and n ==[]:
#             gene_id =  m[0]
#             gene_name = m[0]
#             if gene_id not in id_name_dct:
#                 id_name_dct[gene_id] = gene_name
#     f1.close()
if h2:
    a = 'gene_id=([^;]+);'
    b = 'Name=([^;]+);'
    id_name_dct = get_value(a,b)
else:
    h3 = re.search(".gff", name)
    if h3:
    	a = 'ID=gene-([^;]+);'
    	b = 'GeneID:([^;]+);'
    	id_name_dct = get_value(b,a)
# if h2:
#     for line in f1:
#         if line.startswith("#"):
#             continue
#         line = line.split("\t")[8]

#         m = re.findall(r"gene_id=([^;]+);", line)
#         if m != []:
#             n = re.findall(r"Name=([^;]+);", line)
#             if m != [] and n !=[]:
#                 gene_id = m[0]
#                 gene_name = n[0]
#                 if gene_id not in id_name_dct:
#                     id_name_dct[gene_id] = gene_name
#             if m != [] and n ==[]:
#                 gene_id =  m[0]
#                 gene_name = m[0]
#                 if gene_id not in id_name_dct:
#                     id_name_dct[gene_id] = gene_name
#         if m == []:
#             m = re.findall(r"ID=gene-([^;]+);", line)
#             n = re.findall(r"Name=([^;]+);", line)
#             if m != [] and n !=[]:
#                 gene_id = m[0]
#                 gene_name = n[0]
#                 if gene_id not in id_name_dct:
#                     id_name_dct[gene_id] = gene_name
#     f1.close()
# else:
#     h3 = re.search(".gff", name)
#     if h3 : 
#         for line in f1:
#             if line.startswith("#"):
#                 continue
#             line = line.split("\t")[8]
#             n = re.findall(r"ID=gene-([^;]+);", line)
#             if n != []:
#                 m = re.findall(r"GeneID:([^;]+);", line)
#                 if m != [] and n !=[]:
#                     gene_id = m[0]
#                     gene_name = n[0]
#                     if gene_id not in id_name_dct:
#                         id_name_dct[gene_id] = gene_name
#                 if m != [] and n ==[]:
#                     gene_id =  m[0]
#                     gene_name = m[0]
#                     if gene_id not in id_name_dct:
#                         id_name_dct[gene_id] = gene_name
#         f1.close()


f2 = open(outfile, "w")
f2.write("gene_id" + "\t"+ "gene_name" + "\n")
for k, v in id_name_dct.items():
    line = k + "\t" + v + "\n"
    f2.write(line)
f2.close()