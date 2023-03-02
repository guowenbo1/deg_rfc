#! /nfs1/public2/User/lixiaofei/software/miniconda3/envs/python3/bin/python3
import sys
import os
import argparse
import configparser

# read config file
config = configparser.ConfigParser()
config.read("/nfs3/public2/User/lixiaofei/pipeline/config.cfg")
parser = argparse.ArgumentParser(description="COG Annotation")
parser.add_argument("--query", "-q", type=str, required=True, help="query fasta")
parser.add_argument("--type", "-t", type=str, required=True, choices=['faa', 'fna'], help="type of query file")
parser.add_argument("--visiual_prefix", "-vp", type=str, default=False, help="visual COG with barplot or not")
parser.add_argument("--clear", "-c", action="store_true", default=False, help="clear tmp file or not")
parser.add_argument("--dry_run", "-dr", action="store_true", default=False, help="just output command")
parser.add_argument("--output", "-o", type=str, required=True, help="output file")
args = parser.parse_args()


def file_exists(file):
	if not os.path.exists(file):
		print("%s not exists!" % file)
		exit(1)

outdir = os.path.split(args.visiual_prefix)[0]
if not os.path.exists(outdir):
    os.makedirs(outdir)

file_exists(args.query)

if args.type == 'faa':
	diamond_prog = "blastp"
else:
	diamond_prog = "blastx"

blast_res = "%s.%s.out" % (args.query, diamond_prog)
cmd = """%s %s -q %s -d %s -p 10 -e 1e-5 -k 1 --sensitive -o %s && \\\n""" % (config['Software']['diamond'], diamond_prog, args.query, config['Database']['COG'], blast_res)

cmd += """%s %s -i1 %s -i2 %s -o %s && \\\n""" % (config['Software']['Perl5'], config['COG']['Get_COG_anno'], blast_res, config['COG']['COG_def'], args.output)

if args.visiual_prefix:
	cmd += """%s %s -i %s -p %s.pdf -f %s.png\n""" % (config['Software']['Rscript3.6'], config['COG']['COG_barplot'], args.output, args.visiual_prefix, args.visiual_prefix)

if args.clear:
	cmd += "rm %s" % (blast_res)

if args.dry_run:
	print(cmd)
else:
	os.system(cmd)


