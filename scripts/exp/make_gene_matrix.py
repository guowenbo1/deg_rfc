import os
import pandas as pd
import sys
import re

def get_sample(file):
    sample = os.path.dirname(file)
    sample = os.path.basename(sample)
    return sample


def get_gene_matrix(files, index="FPKM"):
    d = {}
    for file in files:
        sample = get_sample(file)
        df = pd.read_csv(file, sep="\t")
        df = df.ix[:,["Gene ID", index]]
        df = df.groupby("Gene ID").sum()
        s = df[index]
        d[sample] = s
    df = pd.DataFrame(d)
    return df


def get_transcript_matrix(files, index="FPKM"):
    if index == "FPKM":
        num = 4
    elif index == "TPM":
        num = 5
    elif index == "Coverage":
        num = 3
    d = {}
    for file in files:
        f = open(file)
        sample = get_sample(file)
        dic = {}
        for line in f:
            if line.startswith("#"):
                continue
            line_lst = line.split("\t")
            if line_lst[2] != "transcript":
                continue
            data = line_lst[-1]
            # data = data.split(";")
            key = re.findall(r"transcript_id \"(.+?)\";",data)[0]
            if num == 4:
                value = re.findall(r"FPKM \"(.+?)\";",data)[0]
            else:
                value = re.findall(r"TPM \"(.+?)\";",data)[0]
            if key in dic:
                dic[key] += value
            else:
                dic[key] = value
        s = pd.Series(dic)
        d[sample] = s
    df = pd.DataFrame(d)
    return df


gene_files = snakemake.input.gene_abundances
gene_fpkm = get_gene_matrix(index="FPKM", files=gene_files)
gene_tpm = get_gene_matrix(index="TPM", files=gene_files)

gene_fpkm.to_csv(snakemake.output.gene_fpkm, sep="\t")
gene_tpm.to_csv(snakemake.output.gene_tpm, sep="\t")

transcript_files = snakemake.input.gtf
transcript_fpkm = get_transcript_matrix(transcript_files, index="FPKM")
transcript_tpm = get_transcript_matrix(transcript_files, index="TPM")

transcript_fpkm.to_csv(snakemake.output.trans_fpkm, sep="\t")
transcript_tpm.to_csv(snakemake.output.trans_tpm, sep="\t")