#guoyang 2013.11.05
import argparse
from Bio import SeqIO
import gzip
import os.path
import os
# from Bio.Blast.Applications import NcbiblastxCommandline
# from Bio.Blast import NCBIXML

parser = argparse.ArgumentParser(description="blastx source fasta sequences to string protein sequences of specific specie and construct the target protein-protein network with certain score")
parser.add_argument('--species',help="target species code, refer to file: /NJPROJ2/RNA/database/string_ppi/species.v10.txt",required=True)
parser.add_argument('--score',type=int, help="the minimal interaction confidence score",default=400)
parser.add_argument('--fa',help="the query fasta sequences or the total fasta sequences file",required=True)
parser.add_argument('--entries',help="the entries list file that you want to be queried",default=None)
# parser.add_argument('--evalue',help="the E value threshold",default='1e-10')
parser.add_argument('--evalue',help="the E value threshold",default='1e-10')
parser.add_argument('--identity',type=int,help="sequences identity threshold",default=60)
parser.add_argument('--name',default='',help="name of result file")
parser.add_argument('--output',default='',help="dir path to output")
# parser.add_argument('--deg',help="the .DEG.xls file for fold_change and p-value attributes",default=None)
argv=parser.parse_args()
species=argv.species
score=argv.score
fa=argv.fa
entries=argv.entries
e=argv.evalue
identity=argv.identity
name=argv.name
out_name=argv.output
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

print("Extracting sequence for %s" % species)
cur_dir=os.getcwd()
os.chdir('/nfs1/public2/Pipe2/meta/database/string_ppi_v11')
if not os.path.isfile(species+'.protein.sequence.fasta') or os.stat(species+'.protein.sequence.fasta').st_size==0:
    f=gzip.open('protein.sequences.v11.0.fa.gz','rb')
    target_sequence=[]
    for seq_record in SeqIO.parse(f,'fasta'):
        if seq_record.id.startswith(species+'.'):
            target_sequence.append('>'+seq_record.id+'\n')
            target_sequence.append(str(seq_record.seq)+'\n')    
    open(species+'.protein.sequence.fasta','w').writelines(target_sequence)
    f.close()
    print("Extracting sequence for %s done!" % species)
    assert not os.system('/nfs1/public2/Pipe2/meta/software/blast-2.2.26/bin/formatdb -i %s' % (species+'.protein.sequence.fasta'))
target='/nfs1/public2/Pipe2/meta/database/string_ppi_v11/'+species+'.protein.sequence.fasta'
target_dmnd='/nfs1/public2/Pipe2/meta/database/string_ppi_v11/'+species+'.protein.sequence.dmnd'
os.chdir(cur_dir)
#blastx  and  unigene vs reference
# e=float(argv['evalue'])
# blastx_cline = NcbiblastxCommandline(query=query, db=target, evalue=e, max_target_seqs=1, num_threads=10, outfmt=6, out=species+'.blast.txt')
# stdout, stderr = blastx_cline()

if not os.path.exists(target_dmnd):
    os.system("/nfs1/public2/Pipe2/miniconda3/envs/noref-anno/bin/diamond makedb --in %s --db %s" %(target,target_dmnd))
print("Diamond blastp...")
assert not os.system("/nfs1/public2/Pipe2/miniconda3/envs/noref-anno/bin/diamond blastx --query %s --db %s --evalue %s --max-target-seqs 1 --threads 10 --outfmt 6 --out %s.blast.txt" % (query,target_dmnd,e,species))
print("Diamond blastp done!")


name_mapping={}
name_mapping_list=[]
name_mapping_list.append('query id\tsubject id\t%identity\talignment length\tmismatches\tgap opens\tq.start\tq.end\ts.start\ts.end\tevalue\tbit score\n')
for blast_record in open(species+'.blast.txt'):
    if blast_record.strip() != '':
        temp=blast_record.strip().split()
        if float(temp[2]) >= identity:
            name_mapping_list.append(blast_record)
            if temp[1] not in name_mapping:
                name_mapping[temp[1]]=[]
            name_mapping[temp[1]].append(temp[0])
open('entries_name_mapping.txt','w').writelines(name_mapping_list)

#get the ppi for certain species
os.chdir('/nfs1/public2/Pipe2/meta/database/string_ppi_v11')
if not os.path.isfile(species+'.ppi.gz') or os.stat(species+'.ppi.gz').st_size==0:
    print("Extract from All_species for %s" % species)
    assert not os.system('zgrep %s protein.links.detailed.v11.0.txt.gz | gzip > %s.ppi.gz' % ('^'+'"'+species+'\."',species))
    print("Extract from All_species for %s done!" % species)
ppi='/nfs1/public2/Pipe2/meta/database/string_ppi_v11/'+species+'.ppi.gz'
os.chdir(cur_dir)


p=gzip.open(ppi,'rb')
ppi_result=[]
ppi_set=[]
for eachLine in p:
    if eachLine.strip() != '':
        temp=eachLine.strip().split()
        if (temp[0] in name_mapping) and (temp[1] in name_mapping) and (int(temp[-1]) >= score):
            if set([temp[0],temp[1]]) not in ppi_set:
                ppi_set.append(set([temp[0],temp[1]]))
                ppi_result.append(temp[0].replace((species+'.'),'')+'('+','.join(list(set(name_mapping[temp[0]])))+')\t'+temp[1].replace((species+'.'),'')+'('+','.join(list(set(name_mapping[temp[1]])))+')\t'+temp[-1]+'\n')
p.close()

open(out_name,'w').writelines(ppi_result)


