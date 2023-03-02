# expression analysis

rule make_exp_matrix:
    input:
        gene_abundances = expand(assembly_dir + "/ballgrown/{sample}/gene_abundances.tsv", sample=samples),
        gtf = expand(assembly_dir + "/ballgrown/{sample}/transcripts.gtf", sample=samples),
    output:
        gene_fpkm = exp_dir + "/gene_fpkm_all_samples.txt",
        gene_tpm = exp_dir + "/gene_tpm_all_samples.txt",
        trans_fpkm = exp_dir + "/transcript_fpkm_all_samples.txt",
        trans_tpm = exp_dir + "/transcript_tpm_all_samples.txt",
    script:
        "../scripts/exp/make_gene_matrix.py"


rule get_sample_gtf:
    input:
        gtf = expand(assembly_dir + "/ballgrown/{sample}/transcripts.gtf", sample=samples),
    output:
        sample_gtf = assembly_dir + "/ballgrown/sample_lst.txt",
    run:
        with open(output.sample_gtf, "w") as f:
            for file in input.gtf:
                sample = os.path.dirname(file)
                sample = os.path.basename(os.path.dirname(file))
                line = sample + "\t" + file + "\n"
                f.write(line)


rule make_count_matrix:
    input:
        sample_gtf = assembly_dir + "/ballgrown/sample_lst.txt",
    output:
        gene_count = exp_dir + "/gene_count_matrix.txt",
        trans_count = exp_dir + "/transcript_count_matrix.txt",
    shell:
        "/nfs1/public2/User/wzc/miniconda3/envs/rnaseq/bin/prepDE.py -i {input.sample_gtf} -g {output.gene_count} -t {output.trans_count} &&"
        "sed -i 's/,/\t/g' {output.gene_count} && sed -i 's/,/\t/g' {output.trans_count}"


rule get_exp_name_anno:
    input:
        gene_count = exp_dir + "/gene_count_matrix.txt",
        gene_fpkm = exp_dir + "/gene_fpkm_all_samples.txt",
        gene_tpm = exp_dir + "/gene_tpm_all_samples.txt",
        anno = annotation_dir + "/gene_id2name.txt",
    output:
        count_name_anno = exp_dir + "/gene_count_matrix_anno.txt",
        fpkm_name_anno = exp_dir + "/gene_fpkm_all_samples_anno.txt",
        tpm_name_anno = exp_dir + "/gene_tpm_all_samples_anno.txt",
    params:
        Rscript = config["software"]["Rscript"]
    shell:
        "Rscript scripts/exp/exp_name_anno.R {input.gene_count} {input.gene_fpkm} {input.gene_tpm} {input.anno} \
        {output.count_name_anno} {output.fpkm_name_anno} {output.tpm_name_anno}"

rule exp_boxplot: 
    input:
        matrix = exp_dir + "/gene_{index}_all_samples.txt",
        metadata = config["metadata"],
    output:
        fig = expand(exp_dir + "/boxplot/{{index}}_SampleID.{ext}", ext = fig_exts),
        fig1 = expand(exp_dir + "/boxplot/{{index}}_group.{ext}", ext = fig_exts),
    script:
        "../scripts/exp/exp_boxplot.R"


rule exp_density:
    input:
        matrix = exp_dir + "/gene_{index}_all_samples.txt",
        metadata = config["metadata"],
    output:
        fig = expand(exp_dir + "/density/{{index}}_{{group}}.{ext}", ext = fig_exts),
    script:
        "../scripts/exp/exp_density.R"


rule exp_correlation:
    input:
        matrix = exp_dir + "/gene_{index}_all_samples.txt",
        metadata = config["metadata"],
    output:
        corr = exp_dir + "/correlation/{index}_corr.txt",
        fig = expand(exp_dir + "/correlation/{{index}}_corr.{ext}", ext = fig_exts),
        fig1 = expand(exp_dir + "/correlation/{{index}}_corr_cluster.{ext}", ext = fig_exts),
        fig3 = expand(exp_dir + "/correlation/{{index}}_corr2.{ext}", ext = fig_exts),
        fig4 = expand(exp_dir + "/correlation/{{index}}_corr_cluster2.{ext}", ext = fig_exts),
    params:
        height = 8,
        width = 10
    script:
        "../scripts/exp/exp_correlation.R"


rule exp_pca:
    input:
        matrix = exp_dir + "/gene_{index}_all_samples.txt",
        metadata = config["metadata"],
    output:
        pca = exp_dir + "/PCA/{index}_{group}_pca.txt",
        fig = expand(exp_dir + "/PCA/{{index}}_{{group}}_pca.{ext}", ext = fig_exts),
        fig2 = expand(exp_dir + "/PCA/{{index}}_{{group}}_pca_label.{ext}", ext = fig_exts),
        fig3 = expand(exp_dir + "/PCA/{{index}}_{{group}}_pca_ellipse.{ext}", ext = fig_exts),
    script:
        "../scripts/exp/exp_pca.R"


rule exp_pca_apart:
    input:
        matrix = exp_dir + "/gene_{index}_all_samples.txt",
        metadata = config["metadata"],
    params:
        index = "{index}",
        contrast = "{contrast}",
        Rscript = config["software"]["Rscript"],
    output:
        pca = exp_dir + "/PCA/{contrast}/{index}_{group}_pca.txt",
        fig = expand(exp_dir + "/PCA/{{contrast}}/{{index}}_{{group}}_pca.{ext}", ext = fig_exts),
        fig2 = expand(exp_dir + "/PCA/{{contrast}}/{{index}}_{{group}}_pca_label.{ext}", ext = fig_exts),
        fig3 = expand(exp_dir + "/PCA/{{contrast}}/{{index}}_{{group}}_pca_ellipse.{ext}", ext = fig_exts),
    shell:
        "Rscript scripts/exp/exp_pca_apart.R {input.matrix} {input.metadata} {params.index} {params.contrast} {output.pca} {output.fig[0]} {output.fig[1]} {output.fig2[0]} {output.fig2[1]} {output.fig3[0]} {output.fig3[1]}"


rule exp_distribution:
    input:
        matrix = exp_dir + "/gene_fpkm_all_samples.txt",
    output:
        outdir = directory(exp_dir + "/histogram/")
    params:
        Rscript = config["software"]["Rscript"],
    shell:
        "Rscript scripts/exp/exp_histogram.R {input.matrix} {output.outdir}"


rule exp_upset:
    input:
        matrix = exp_dir + "/gene_count_matrix.txt",
    output:
        upset1 = exp_dir  + "/upset/sample_upset.pdf",    
        upset = exp_dir  + "/upset/sample_upset.png",
    shell: 
        "/nfs1/public2/User/wzc/miniconda3/envs/amplicon/bin/Rscript scripts/exp/exp_upset.R -f {input.matrix} -p {output.upset1} -o {output.upset} "