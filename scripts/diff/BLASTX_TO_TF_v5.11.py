#guoyang 2013.11.05
import argparse
from Bio import SeqIO
import gzip
import os.path
import os
# from Bio.Blast.Applications import NcbiblastxCommandline
# from Bio.Blast import NCBIXML

parser = argparse.ArgumentParser(description="blastx source fasta sequences to string protein sequences of specific specie and construct the target protein-protein network with certain score")
parser.add_argument('--fa',help="the query fasta sequences or the total fasta sequences file",required=True)
parser.add_argument('--entries',help="the entries list file that you want to be queried",default=None)
parser.add_argument('--target',help="the entries list file that you want to be queried",default=None)
parser.add_argument('--TF',help="TF list",default=None)
# parser.add_argument('--evalue',help="the E value threshold",default='1e-10')
parser.add_argument('--evalue',help="the E value threshold",default='1e-5')
parser.add_argument('--identity',type=int,help="sequences identity threshold",default=60)
parser.add_argument('--name',default='',help="name of result file")
parser.add_argument('--output',default='',help="dir path to output")
# parser.add_argument('--deg',help="the .DEG.xls file for fold_change and p-value attributes",default=None)
argv=parser.parse_args()
fa=argv.fa
entries=argv.entries
e=argv.evalue
identity=argv.identity
name=argv.name
out_name=argv.output
target = argv.target
TF = argv.TF
out=os.path.dirname(out_name)
if not os.path.exists(out):
    os.mkdir(out)
assert os.path.isdir(out)
os.chdir(out)
if entries != None:
    entry=[each.strip() for each in open(entries) if each.strip() != '']
    query_sequence=[]
    for seq_record in SeqIO.parse(fa,'fasta'):
        if seq_record.id.split("_i")[0] in entry:
            query_sequence.append('>'+seq_record.id+'\n')
            query_sequence.append(str(seq_record.seq)+'\n')
    query='query_sequences.fasta'
    open(query,'w').writelines(query_sequence)
else:
    query=fa


target_dmnd=os.path.splitext(target)[0] + ".dmnd"
print(target_dmnd)
cur_dir=os.getcwd()
#blastx  and  unigene vs reference
# e=float(argv['evalue'])
# blastx_cline = NcbiblastxCommandline(query=query, db=target, evalue=e, max_target_seqs=1, num_threads=10, outfmt=6, out=species+'.blast.txt')
# stdout, stderr = blastx_cline()

if not os.path.exists(target_dmnd):
    os.system("/nfs1/public2/Pipe2/miniconda3/envs/noref-anno/bin/diamond makedb --in %s --db %s" %(target,target_dmnd))
print("Diamond blastp...")
assert not os.system("/nfs1/public2/Pipe2/miniconda3/envs/noref-anno/bin/diamond blastx --query %s --db %s --evalue %s --max-target-seqs 1 --threads 10 --outfmt 6 --out blast.txt" % (query,target_dmnd,e))
print("Diamond blastp done!")

name_mapping={}
with open(TF, 'r') as in_f_1:
    next(in_f_1)
    for i in in_f_1:
        name_mapping[i.strip().split("\t")[0]] = i.strip().split("\t")[-1]


blast_dct = {}
with open("blast.txt", "r") as in_f_2:
    for i in in_f_2:
        blast_dct[i.strip().split("\t")[0]] = i.strip().split("\t")[1]

res_lst = []
for k, v in blast_dct.items():
    for m, n in name_mapping.items():
        if v == m:
            res_lst.append("%s\t%s\n" % (k, n))

open(out_name,'w').writelines(res_lst)

