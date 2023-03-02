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
        diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls"
    params:
        logfc = config["diff_cutoff"]["logfc"],
        pvalue = config["diff_cutoff"]["pvalue"],
        fdr = config["diff_cutoff"]["fdr"]
    script:
        "../scripts/diff/add_fpkm_tpm_annotation.R"


rule add_all_annotation_diff:
    input:
        diff_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        anno = annotation_dir + "/gene2function.xls"
    output:
        diff_anno_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_anno.xls",
    shell:
        "csvtk join -t -l  {input} > {output}"

rule diff_plot:
    input:
        diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls"
    output:
        fig1 = diff_dir + "/{level}/{method}/{contrast}/{contrast}_Volcano.{ext}",
        fig2 = diff_dir + "/{level}/{method}/{contrast}/{contrast}_MA.{ext}",
    params:
        pvalue = config["diff_cutoff"]["pvalue"]
    shell:
        "Rscript scripts/diff/diff_plot.R {wildcards.method} {wildcards.contrast} {input} {params.pvalue} {output.fig1} {output.fig2}"


rule diff_heatmap:
    input:
        fpkm = exp_dir + "/{level}_fpkm_all_samples.txt",
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        metadata = config["metadata"],
        diff_lst = config["diff_lst"]
    output:
        fig = diff_dir + "/{level}/{method}/{contrast}/{contrast}_heatmap.{ext}",
    script:
        "../scripts/diff/diff_heatmap.R"


rule go_enrich:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        go_anno = annotation_dir + "/go/go.txt",
    output:
        res = diff_dir + "/{level}/{method}/{contrast}/go_enrich/go_enrich_{label}.xls",
        barplot = diff_dir + "/{{level}}/{{method}}/{{contrast}}/go_enrich/barplot_{{label}}.{ext}", 
        dotplot = diff_dir + "/{{level}}/{{method}}/{{contrast}}/go_enrich/dotplot_{{label}}.{ext}",
        DAGpdf_BP_pdf = diff_dir + "/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_BP.pdf",
        DAGpng_BP = diff_dir + "/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_BP.png",
        DAGpdf_CC_pdf = diff_dir + "/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_CC.pdf",
        DAGpng_CC = diff_dir + "/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_CC.png",
        DAGpdf_MF_pdf = diff_dir + "/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_MF.pdf",
        DAGpng_MF = diff_dir + "/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_MF.png" 

    shell:  """

            Rscript scripts/diff/go_enrich.R {wildcards.label} {wildcards.contrast} {input.deg} {input.go_anno} \
            {output.res} {output.barplot} {output.dotplot} {diff_dir}/{level}/{method}/{contrast}/go_enrich/DAG/{label};

            convert {output.DAGpdf_BP_pdf} {output.DAGpng_BP};
            convert {output.DAGpdf_CC_pdf} {output.DAGpng_CC};
            convert {output.DAGpdf_MF_pdf} {output.DAGpng_MF}

            """



rule kegg_enrich:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        kegg_anno = annotation_dir + "/kegg/kegg.txt",
    output:
        res = diff_dir + "/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_{label}.xls",
        barplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/kegg_enrich/barplot_{{label}}.{ext}", ext=fig_exts),
        dotplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/kegg_enrich/dotplot_{{label}}.{ext}", ext=fig_exts)
    script:
        "../scripts/diff/kegg_enrich.R"
