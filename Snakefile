# ref mRNA-seq analysis
import os
import pandas as pd

# load configfile
from snakemake.io import expand

configfile: "config.yaml"

# global vals
metadata = pd.read_csv(config["metadata"], sep="\t").set_index("SampleID", drop=False)

samples = metadata["SampleID"]
samples_num = len(samples)
species = config["species"]

groups = config["groups"]
diff_levels = config["diff_levels"]
diff_methods = config["diff_method"]

# 重构contrast的获取，方便修改

diff_lst = config["diff_lst"]
diff_groups = []
contrasts = []
contrasts_groups = {}
with open(diff_lst, "r") as f:
    for line in f:
        line_list = line.split()
        contrasts.append(line_list[0])
        diff_groups.append(line_list[1])
        contrasts_groups[line_list[0]] = [line_list[2], line_list[3]]

fig_exts = ["pdf", "png"]
deg_labels = ["all", "up", "down"]
types = ["all_weblink"]

# dirs
raw_data_dir = config["dirs"][0]
clean_data_dir = config["dirs"][1]
align_dir = config["dirs"][2]
assembly_dir = config["dirs"][3]
exp_dir = config["dirs"][4]
new_trans_dir = config["dirs"][5]
diff_dir = config["dirs"][6]
as_dir = config["dirs"][7]
annotation_dir = config["dirs"][8]
seqkit_stats_dir = config["dirs"][9]


# wildcard_constraints
wildcard_constraints:
    sample = "[a-zA-Z][a-zA-Z0-9_]*",
    group = "[a-zA-Z][a-zA-Z0-9_]*",

# rule target
rule all:
    input:
        # clean_data
        # expand(seqkit_stats_dir + "/{sample}.xls", sample=samples),
        # clean_data_dir + "/fastp/fastp_stats.xls",
        # align_dir + "/hisat2/hisat2_stats.xls",
        # clean_data_dir + "/filter_stats.xls",
        # expand(clean_data_dir + "/fastp/{sample}/{sample}.GC.png", sample=samples),
        # expand(clean_data_dir + "/fastp/{sample}/{sample}.Error.png", sample=samples),
        # # Map/QC
        # expand(align_dir + "/QC/read_distribution/{sample}.read_distribution.txt", sample=samples),
        # expand(align_dir + "/QC/read_distribution/{sample}.read_distribution_stat.txt", sample=samples),
        # expand(align_dir + "/QC/read_distribution/{sample}.read_distribution_plot.pdf", sample=samples),
        # #expand(align_dir + "/QC/junctionSaturation/{sample}.junctionSaturation_plot.pdf", sample=samples),
        expand(align_dir + "/QC/geneBody_coverage/{sample}.geneBodyCoverage.curves.pdf", sample=samples),
        # expand(align_dir + "/QC/read_duplication/{sample}.DupRate_plot.pdf", sample=samples),
        # expand(expand(align_dir + "/QC/RPKM_saturation/{{sample}}_saturation_curve.{ext}", ext=fig_exts), sample=samples),
        # expand(assembly_dir + "/ballgrown/{sample}/transcripts.gtf", sample=samples),

        # Assembly
        "Assembly/comparison/",

        # AS
        new_trans_dir + "/new_trans.gtf",

        # Annotation
        annotation_dir + "/genes.fa",
        annotation_dir + "/eggnog/eggnog.txt",
        annotation_dir + "/go/go.txt",
        annotation_dir + "/kegg/kegg.txt",
        annotation_dir + "/COG/COG_annotations.xls",

        # Expression
        exp_dir + "/gene_fpkm_all_samples.txt",
        exp_dir + "/gene_tpm_all_samples.txt",
        exp_dir + "/transcript_fpkm_all_samples.txt",
        exp_dir + "/transcript_tpm_all_samples.txt",
        exp_dir + "/gene_count_matrix.txt",
        exp_dir + "/transcript_count_matrix.txt",
        exp_dir + "/gene_count_matrix_anno.txt",
        exp_dir + "/gene_fpkm_all_samples_anno.txt",
        exp_dir + "/gene_tpm_all_samples_anno.txt",
        expand(expand(exp_dir + "/boxplot/{{index}}_{{group}}.{ext}", ext=fig_exts), index=["fpkm", "tpm"], group=config["groups"]+["SampleID"]),
        expand(expand(exp_dir + "/density/{{index}}_{{group}}.{ext}", ext=fig_exts), index=["fpkm", "tpm"], group=config["groups"]+["SampleID"]),
        expand(exp_dir + "/correlation/{index}_corr.txt", index=["fpkm", "tpm"]),
        expand(exp_dir + "/PCA/{index}_{group}_pca.txt", index=["fpkm", "tpm"], group=config["groups"]),
        expand(expand(exp_dir + "/PCA/{{index}}_{{group}}_pca.{ext}", ext=fig_exts), index=["fpkm", "tpm"], group=config["groups"]),
        expand(expand(exp_dir + "/PCA/{{index}}_{{group}}_pca_label.{ext}", ext=fig_exts), index=["fpkm", "tpm"], group=config["groups"]),
        expand(expand(exp_dir + "/PCA/{{index}}_{{group}}_pca_ellipse.{ext}", ext=fig_exts), index=["fpkm", "tpm"], group=config["groups"]),
        # exp_apart
        expand(expand(exp_dir + "/PCA/{{contrast}}/{{index}}_{{group}}_pca.{ext}", ext=fig_exts), index=["fpkm", "tpm"], group=config["groups"], contrast=contrasts),
        expand(expand(exp_dir + "/PCA/{{contrast}}/{{index}}_{{group}}_pca_label.{ext}", ext=fig_exts), index=["fpkm", "tpm"], group=config["groups"], contrast=contrasts),
        expand(expand(exp_dir + "/PCA/{{contrast}}/{{index}}_{{group}}_pca_ellipse.{ext}", ext=fig_exts), index=["fpkm", "tpm"], group=config["groups"], contrast=contrasts),
        expand(exp_dir + "/histogram/"),
        exp_dir+ "/upset/sample_upset.png",

        # Diff
        expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_anno.xls", level=diff_levels, contrast=contrasts, method=diff_methods),
        expand(diff_dir + "/{level}/{method}/DEGs_venn_or_flower.pdf",level=diff_levels, method=diff_methods),
        expand(diff_dir + '/{level}/{method}/All_DEGs/',level=diff_levels, method=diff_methods),          
        expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_Volcano.{ext}", level=diff_levels, contrast=contrasts, method=diff_methods, ext=fig_exts),
        expand(diff_dir + "/{level}/{method}/all_diff_plot.png",level=diff_levels, method=diff_methods), 
        expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_MA.{ext}", level=diff_levels, contrast=contrasts, method=diff_methods, ext=fig_exts),
        expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_heatmap.{ext}", level=diff_levels, contrast=contrasts, method=diff_methods, ext=fig_exts),
        expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_heatmap_details.pdf", level=diff_levels, contrast=contrasts, method=diff_methods, ext=fig_exts),
        expand(diff_dir + "/{level}/{method}/{contrast}/go_enrich/go_enrich_{label}.xls", level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),
        expand(expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/go_enrich/barplot_{{label}}.{ext}", ext=fig_exts),level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),
        expand(expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/go_enrich/dotplot_{{label}}.{ext}", ext=fig_exts),level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),
        expand(diff_dir + "/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_BP.pdf",level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),
        expand(diff_dir + "/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_CC.pdf",level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),
        expand(diff_dir + "/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_MF.pdf",level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),
        expand(diff_dir + "/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_{label}.xls", level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),
        expand(expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/kegg_enrich/barplot_{{label}}.{ext}", ext=fig_exts),level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),
        expand(expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/kegg_enrich/dotplot_{{label}}.{ext}", ext=fig_exts),level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),
        expand(diff_dir + '/{level}/{method}/{contrast}/gsea/',level=diff_levels, contrast=contrasts, method=diff_methods),
        expand(diff_dir + '/{level}/{method}/{contrast}/cog/',level=diff_levels, contrast=contrasts, method=diff_methods),
        # cricle
        expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_circle.png", level=diff_levels, contrast=contrasts, method=diff_methods),

        # # PPI
        #expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_PPI.xls", level=diff_levels, contrast=contrasts, method=diff_methods),
        # # TF
        #expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_TF.xls", level=diff_levels, contrast=contrasts, method=diff_methods, label=deg_labels),

        # # Report
        #"mRNA_report.html",


        # Deprecated
        # expand(diff_dir + "/{level}/{method}/{contrast}/kegg_enrich/maps/",level=diff_levels, contrast=contrasts, method=diff_methods),
        #expand(diff_dir + "/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_{type}.xls",level=diff_levels, contrast=contrasts, method=diff_methods, type=types),        
        # expand(as_dir + "/{contrast}/results", contrast=contrasts),
        
# includes
include: "rules/trim.py"
include: "rules/align.py"
include: "rules/assembly.py"
include: "rules/new_trans.py"
include: "rules/annotation.py"
include: "rules/exp.py"
include: "rules/diff.py"
include: "rules/report.py"
include: "rules/rmats.py"
