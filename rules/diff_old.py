# diff analysis with trinity script

rule get_samples_file:
    input:
        config["metadata"]
    output:
        diff_dir + "/samples_file.txt"
    params:
        group = groups
    run:
        df = pd.read_csv(input[0], sep="\t")
        i = 0 
        for group in params.group: 
            if i > 0: 
                temp = data
                data = df.loc[:, [group, "SampleID"]] 
                data.columns = ["group", "SampleID"] 
                data = pd.concat([data, temp])  
            else: 
                data = df.loc[:, [group, "SampleID"]] 
                data.columns = ["group", "SampleID"] 
                i += 1 
        data.to_csv(output[0], index=False, header=False, sep='\t')


rule get_contrasts:
    input:
        config["diff_lst"]
    output:
        diff_dir + "/contrasts.txt"
    run:
        with open(output[0], "w") as f:
            with open(input[0]) as f_in:
                for line in f_in:
                    line = line.strip("\n")
                    line_lst = line.split("\t")
                    line = line_lst[2] + "\t" + line_lst[3] + "\n"
                    f.write(line)


rule diff:
    input:
        matrix = exp_dir + "/{level}_count_matrix.txt",
        samples_file = diff_dir + "/samples_file.txt",
        contrasts = diff_dir + "/contrasts.txt",
    output:
        # directory(diff_dir + "/{level}"),
        expand(diff_dir + "/{{level}}/{contrast}.{{method}}.DE_results", contrast=contrasts),
        expand(diff_dir + "/{{level}}/{contrast}.{{method}}.count_matrix", contrast=contrasts),
    params:
        # method = config["diff_method"],
        dispersion = 0.1
    run:
        method = wildcards.method
        outdir = os.path.dirname(output[0])
        cmd = "run_DE_analysis.pl --matrix {input.matrix} --samples_file {input.samples_file} --contrasts {input.contrasts} --output {outdir} "
        if method == "DESeq2":
            cmd += " --method DESeq2 "
        elif method == "edgeR":
            cmd += " --method edgeR --dispersion {params.dispersion} "
        cmd += "&& sed -i 's/sampleA/{wildcards.level}\tsampleA/' {outdir}/*{method}.DE_results"
        cmd += "&& sed -i '1s/^/{wildcards.level}\t/' {outdir}/*{method}.count_matrix"
        shell(cmd)
        prefix = os.path.basename(input.matrix)
        for file in os.listdir(outdir):
            file_new = file.replace(prefix+".", "")
            os.rename(outdir + "/" + file, outdir + "/" + file_new)


rule add_fpkm_tpm_annotation:
    input:
        de_res = diff_dir + "/{level}/{contrast}.{method}.DE_results",
        fpkm = exp_dir + "/{level}_fpkm_all_samples.txt",
        tpm = exp_dir + "/{level}_tpm_all_samples.txt",
        gene_id2name = annotation_dir + "/gene_id2name.txt"
    output:
        diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls"
    params:
        logfc = config["diff_cutoff"]["logfc"],
        pvalue = config["diff_cutoff"]["pvalue"],
        fdr = config["diff_cutoff"]["fdr"]
    script:
        "../scripts/diff/add_fpkm_tpm_annotation.R"


rule add_all_annotation_diff:
    input:
        diff_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        anno = annotation_dir + "/gene2function.xls",
        anno_temp = annotation_dir + "/gene2function_temp.xls"
    output:
        diff_anno_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all_anno.xls",
        diff_anno_tab_temp = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all_anno_temp.xls"
    shell:
        "csvtk join -t -l -k {input.diff_tab} {input.anno} > {output.diff_anno_tab} && csvtk join -t -l {input.diff_tab} {input.anno_temp} > {output.diff_anno_tab_temp}"


rule break_diffanno:
    input:
        diff_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        diff_anno_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all_anno.xls",
        diff_anno_tab_temp = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all_anno_temp.xls",
    output:
        diff_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        diffanno = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_anno.xls",
        diffanno_temp = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_anno_temp.xls",
    shell:
        "Rscript scripts/diff/diff_anno.R {input.diff_tab} {output.diff_tab} && Rscript scripts/diff/diff_anno.R {input.diff_anno_tab} {output.diffanno} && Rscript scripts/diff/diff_anno.R {input.diff_anno_tab_temp} {output.diffanno_temp}"
        # "Rscript scripts/diff/diff_anno.R {input.diff_anno_tab} {output.diffanno}"
#各组差异基因venn图或upset图
rule DEGs_list:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls"
    output:
        DEGs = diff_dir + "/{level}/{method}/DEGs_list/{contrast}_DEGs_list.xls"
    params:
        contrast = "{contrast}"
    shell:
        "Rscript scripts/diff/DEG_list.R  {input.deg} {params.contrast} {output.DEGs} "

rule DEGs_all_list:
    input:
        expand(diff_dir + "/{level}/{method}/DEGs_list/{contrast}_DEGs_list.xls", level=diff_levels, contrast=contrasts, method=diff_methods)
    output:
        diff_dir + "/{level}/{method}/DEGs_list/all_DEGs_list.xls"
    run:
        outdir = os.path.dirname(input[0])
        shell("paste  {outdir}/* > {output} ")

rule diff_veen:
    input:
        diff_dir + "/{level}/{method}/DEGs_list/all_DEGs_list.xls"
    output:
        fig1 = diff_dir + "/{level}/{method}/DEGs_venn_or_flower.pdf",
        fig2 = diff_dir + "/{level}/{method}/DEGs_venn_or_flower.png",
    shell:
        "/nfs1/public2/User/lixiaofei/software/miniconda3/envs/R/bin/Rscript scripts/diff/diff_venn.R -f {input} -p {output.fig1} -o {output.fig2}"

#差异基因火山图，MA图
rule diff_plot:
    input:
        diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls"
    output:
        fig1 = diff_dir + "/{level}/{method}/{contrast}/{contrast}_Volcano.{ext}",
        fig2 = diff_dir + "/{level}/{method}/{contrast}/{contrast}_MA.{ext}",
    params:
        pvalue = config["diff_cutoff"]["pvalue"]
    shell:
        "Rscript scripts/diff/diff_plot.R {wildcards.method} {wildcards.contrast} {input} {params.pvalue} {output.fig1} {output.fig2}"

#差异基因热图
rule diff_heatmap:
    input:
        fpkm = exp_dir + "/{level}_fpkm_all_samples.txt",
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        metadata = config["metadata"],
        diff_lst = config["diff_lst"]
    output:
        fig = diff_dir + "/{level}/{method}/{contrast}/{contrast}_heatmap.{ext}",
    script:
        "../scripts/diff/diff_heatmap.R"


rule go_enrich:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        go_anno = annotation_dir + "/go/go.txt",
    output:
        res = diff_dir + "/{level}/{method}/{contrast}/go_enrich/go_enrich_{label}.xls",
        barplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/go_enrich/barplot_{{label}}.{ext}", ext=fig_exts),
        dotplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/go_enrich/dotplot_{{label}}.{ext}", ext=fig_exts)
    script:
        "../scripts/diff/go_enrich.R" 

rule go_dag:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        go_anno = annotation_dir + "/go/go.txt",
    output:
        BPpdf = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_BP.pdf",
        BPpng = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_BP.png",
        CCpdf = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_CC.pdf",
        CCpng = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_CC.png",
        MFpdf = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_MF.pdf",
        MFpng = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_MF.png"
    # shell:
    #      "Rscript scripts/diff/go_dag.R {wildcards.label} {wildcards.contrast} {input.deg} {input.go_anno} {diff_dir}/{level}/{method}/{contrast}/go_enrich/DAG/{label} &&"
    #      "convert {output.BPpdf} {output.BPpng} &&"
    #      "convert {output.CCpdf} {output.CCpng} &&"
    #      "convert {output.MFpdf} {output.MFpng} "
    run:
        outdir = os.path.dirname(output[0])
        cmd = "/nfs1/public2/User/wzc/miniconda3/envs/rnaseq/bin/Rscript scripts/diff/go_dag.R {wildcards.label} {wildcards.contrast} {input.deg} {input.go_anno} {outdir} &&"
        cmd += "convert {output.BPpdf} {output.BPpng} &&"
        cmd += "convert {output.CCpdf} {output.CCpng} &&"
        cmd += "convert {output.MFpdf} {output.MFpng} "        
        shell(cmd)

rule kegg_enrich:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        kegg_anno = annotation_dir + "/kegg/kegg.txt",
    output:
        res = diff_dir + "/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_{label}.xls",
        barplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/kegg_enrich/barplot_{{label}}.{ext}",ext=fig_exts),
        dotplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/kegg_enrich/dotplot_{{label}}.{ext}",ext=fig_exts)
    script:
        "../scripts/diff/kegg_enrich.R"    


# new_add
rule PPI:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        fasta = config['ref']['protein']
    params:
        species_code = config['species_code'],
        contrast = '{contrast}'
    output:
        temp = temp(diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_gene_lst.xls"),
        PPI_res = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_PPI.xls"
    run:
        cwd = os.getcwd()
        cmd = """awk -F '\\t' '{{if($13 == "Decreased" || $13 == "Increased") print $12}}' {input.deg} > {output.temp} && """
        cmd += "/nfs1/public2/Pipe2/miniconda3/envs/noref-ppi/bin/python scripts/diff/BLASTX_TO_PPI_v5.11.py --entries %s/{output.temp} --species {params.species_code} --fa {input.fasta} --name %s/{params.contrast} --output %s/{output.PPI_res}" % (cwd, cwd, cwd)
        shell(cmd)

rule TF:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        fasta = config['ref']['protein'],
        target = "/nfs2/public2/User/lixiaofei/projects/mRNA/20210720-HT20210621236-chenguangxia/genome_file/Vvi_pep.fas",
        TF_list = "/nfs2/public2/User/lixiaofei/projects/mRNA/20210720-HT20210621236-chenguangxia/genome_file/Vvi_TF_list.txt"
    params:
        contrast = '{contrast}'
    output:
        temp = temp(diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_gene_lst.xls"),
        TF_res = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_TF.xls"
    run:
        cwd = os.getcwd()
        cmd = """awk -F '\\t' '{{if($13 == "Decreased" || $13 == "Increased") print $12}}' {input.deg} > {output.temp} && """
        cmd += "/nfs1/public2/Pipe2/miniconda3/envs/noref-ppi/bin/python scripts/diff/BLASTX_TO_TF_v5.11.py --entries %s/{output.temp} --fa {input.fasta} --TF {input.TF_list} --name %s/{params.contrast} --target {input.target} --output %s/{output.TF_res}" % (cwd, cwd, cwd)
        shell(cmd)


rule ggbio_circle:
    input:
        fpkm = exp_dir + "/gene_fpkm_all_samples_anno.txt",
        metadata = config["metadata"],
        chrom_info = config["ref"]["chrom_info"],
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
    params:
        contrast = "{contrast}",  
        ggbio_parameter = config['ggbio_parameter'],
    output:
        circle_pdf = diff_dir + "/{level}/{method}/{contrast}/{contrast}_circle.pdf",
        circle_png = diff_dir + "/{level}/{method}/{contrast}/{contrast}_circle.png",
    shell:
        "/nfs1/public2/User/axl/miniconda3/bin/Rscript scripts/diff/ggbio.R {input.fpkm} {input.metadata} {input.chrom_info} {params.contrast} {input.deg} {params.ggbio_parameter} {output.circle_pdf} {output.circle_png} "