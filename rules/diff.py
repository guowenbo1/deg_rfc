# diff analysis with trinity script

rule get_samples_file:
    input:
        config["metadata"],
    output:
        diff_dir + "/samples_file.txt",
    params:
        group = groups,
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
        config["diff_lst"],
    output:
        diff_dir + "/contrasts.txt",
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
        expand(diff_dir + "/{{level}}/{contrast}.{{method}}.DE_results", contrast=contrasts),
        expand(diff_dir + "/{{level}}/{contrast}.{{method}}.count_matrix", contrast=contrasts),
    params:
        run_DE = config["software"]["run_DE"],
        method = config["diff_method"],
        dispersion = 0.1
    run:
        outdir = os.path.dirname(output[0])
        cmd = "{params.run_DE} --matrix {input.matrix} --samples_file {input.samples_file} --contrasts {input.contrasts} --output {outdir} "
        diff_method = str(params.method)
        if diff_method == "DESeq2":
            cmd += " --method DESeq2 "
        elif diff_method == "edgeR":
            cmd += " --method edgeR --dispersion {params.dispersion} "
        cmd += "&& sed -i 's/sampleA/{wildcards.level}\tsampleA/' {outdir}/*" + diff_method + ".DE_results"
        cmd += "&& sed -i '1s/^/{wildcards.level}\t/' {outdir}/*" + diff_method + ".count_matrix"
        print(cmd)
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
        gene_id2name = annotation_dir + "/gene_id2name.txt",
    output:
        diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
    params:
        type = config["diff_cutoff"]["type"],
        value = config["diff_cutoff"]["value"],
        logfc = config["diff_cutoff"]["logfc"],
    script:
        "../scripts/diff/add_fpkm_tpm_annotation.R"


rule add_all_annotation_diff:
    input:
        diff_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        anno = annotation_dir + "/gene2function.xls",
        anno_temp = annotation_dir + "/gene2function_temp.xls",
    output:
        diff_anno_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all_anno.xls",
        diff_anno_tab_temp = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all_anno_temp.xls",
    params:
        csvtk = config["software"]["csvtk"]
    shell:
        "{params.csvtk} join -t -l -k {input.diff_tab} {input.anno} > {output.diff_anno_tab} && {params.csvtk} join -t -l {input.diff_tab} {input.anno_temp} > {output.diff_anno_tab_temp}"


rule break_diffanno:
    input:
        diff_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        diff_anno_tab = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all_anno.xls",
        diff_anno_tab_temp = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all_anno_temp.xls",
    output:
        diff = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        diffanno = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_anno.xls",
        diffanno_temp = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_anno_temp.xls",
    params:
        Rscript = config["software"]["Rscript"]
    shell:
        "Rscript scripts/diff/diff_anno.R {input.diff_anno_tab} {output.diffanno} && Rscript scripts/diff/diff_anno.R {input.diff_tab} {output.diff} && Rscript scripts/diff/diff_anno.R {input.diff_anno_tab_temp} {output.diffanno_temp}"


#各组差异基因venn图或upset图
rule DEGs_list:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
    output:
        DEGs = diff_dir + "/{level}/{method}/DEGs_list/{contrast}_DEGs_list.xls",
    params:
        contrast = "{contrast}",
        Rscript = config["software"]["Rscript"]
    shell:
        "Rscript scripts/diff/DEG_list.R {input.deg} {params.contrast} {output.DEGs}"


rule DEGs_all_list:
    input:
        expand(diff_dir + "/{level}/{method}/DEGs_list/{contrast}_DEGs_list.xls", level=diff_levels, contrast=contrasts, method=diff_methods),
    output:
        diff_dir + "/{level}/{method}/DEGs_list/all_DEGs_list.xls",
    run:
        # dirs = os.path.dirname(input[0])
        shell("paste {input} > {output}")

rule diff_barplot:
    input:
        expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls", level=diff_levels, contrast=contrasts, method=diff_methods)
    output:
        out = directory(diff_dir + '/{level}/{method}/All_DEGs/'),
    shell:
        "/nfs1/public2/User/axl/miniconda3/bin/Rscript scripts/diff/stat.DEG.num.plot.R -i {input} --output {output.out} --name All.DEGs.stat --name All.DEGs.stat "

rule diff_veen:
    input:
        diff_dir + "/{level}/{method}/DEGs_list/all_DEGs_list.xls",
    output:
        fig1 = diff_dir + "/{level}/{method}/DEGs_venn_or_flower.pdf",
        fig2 = diff_dir + "/{level}/{method}/DEGs_venn_or_flower.png",
    shell:
        "/nfs1/public2/User/lixiaofei/software/miniconda3/envs/R/bin/Rscript scripts/diff/diff_venn.R -f {input} -p {output.fig1} -o {output.fig2} &> /dev/null"


#差异基因火山图，MA图
rule diff_plot:
    input:
        diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
    output:
        fig1 = diff_dir + "/{level}/{method}/{contrast}/{contrast}_Volcano.{ext}",
        fig2 = diff_dir + "/{level}/{method}/{contrast}/{contrast}_MA.{ext}",
    params:
        type = config["diff_cutoff"]["type"],
        value = config["diff_cutoff"]["value"],
        logfc = config["diff_cutoff"]["logfc"],
        Rscript = config["software"]["Rscript"],
    shell:
        "Rscript scripts/diff/diff_plot.R {wildcards.method} {wildcards.contrast} {input} {params.type} {params.logfc} {params.value} {output.fig1} {output.fig2}"

rule all_diff_plot:
    input:
        expand(diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",level=diff_levels, contrast=contrasts, method=diff_methods),
    output:
        fig1 = diff_dir + "/{level}/{method}/all_diff_plot.pdf",
        fig2 = diff_dir + "/{level}/{method}/all_diff_plot.png",
    params:
        type = config["diff_cutoff"]["type"],
        value = config["diff_cutoff"]["value"],
        logfc = config["diff_cutoff"]["logfc"],
        Rscript = config["software"]["Rscript"],
    shell:
        "/nfs1/public2/User/axl/miniconda3/bin/Rscript scripts/diff/all_diff_plot.R  -i {input} -p {output.fig1} -o {output.fig2}"
        
#差异基因热图
rule diff_heatmap:
    input:
        fpkm = exp_dir + "/{level}_fpkm_all_samples.txt",
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        metadata = config["metadata"],
        diff_lst = config["diff_lst"]
    output:
        fig2 = diff_dir + "/{level}/{method}/{contrast}/{contrast}_heatmap_details.{ext}", 
    script:
        "../scripts/diff/diff_heatmap.R"
rule diff_heatmap2:
    input:
        fpkm = exp_dir + "/{level}_fpkm_all_samples.txt",
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        metadata = config["metadata"],
        diff_lst = config["diff_lst"]
    output:
        fig1 = diff_dir + "/{level}/{method}/{contrast}/{contrast}_heatmap.{ext}", 
    script:
        "../scripts/diff/diff_heatmap2.R"



rule go_enrich:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        go_anno = annotation_dir + "/go/go.txt",
    output:
        res = diff_dir + "/{level}/{method}/{contrast}/go_enrich/go_enrich_{label}.xls",
        barplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/go_enrich/barplot_{{label}}.{ext}", ext=fig_exts),
        dotplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/go_enrich/dotplot_{{label}}.{ext}", ext=fig_exts),
    script:
        "../scripts/diff/go_enrich.R" 

rule go_dag:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        go_anno = annotation_dir + "/go/go.txt",
    output:
        BPpdf = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_BP.pdf",
        BPpng = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_BP.png",
        CCpdf = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_CC.pdf",
        CCpng = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_CC.png",
        MFpdf = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_MF.pdf",
        MFpng = diff_dir+"/{level}/{method}/{contrast}/go_enrich/DAG/{label}/DAGplot_MF.png",
    run:
        outdir = os.path.dirname(output[0])
        cmd = "/nfs1/public2/User/wzc/miniconda3/envs/rnaseq/bin/Rscript scripts/diff/go_dag.R {wildcards.label} {wildcards.contrast} {input.deg} {input.go_anno} {outdir} &&"
        cmd += "convert {output.BPpdf} {output.BPpng} &&"
        cmd += "convert {output.CCpdf} {output.CCpng} &&"
        cmd += "convert {output.MFpdf} {output.MFpng} "        
        shell(cmd)


rule kegg_enrich:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
        kegg_anno = annotation_dir + "/kegg/kegg.txt",
    output:
        res = diff_dir + "/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_{label}.xls",
        barplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/kegg_enrich/barplot_{{label}}.{ext}",ext=fig_exts),
        dotplot = expand(diff_dir + "/{{level}}/{{method}}/{{contrast}}/kegg_enrich/dotplot_{{label}}.{ext}",ext=fig_exts)
    script:
        "../scripts/diff/kegg_enrich.R"
        
rule gsea:
    input:
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
        kegg_anno = annotation_dir + "/kegg/kegg.txt",
        go_anno = annotation_dir + "/go/go.txt",
    output:
        resdir = directory(diff_dir + '/{level}/{method}/{contrast}/gsea/'),
        # out = expand(diff_dir + '/{level}/{method}/{contrast}/gsea/go.xls',level=diff_levels, contrast=contrasts, method=diff_methods),
    shell:
        "/nfs1/public2/User/axl/miniconda3/bin/Rscript scripts/diff/gsea.R  {input.go_anno} {input.kegg_anno}  {input.deg} {output.resdir}go.xls {output.resdir}kegg.xls {output.resdir}"

# TODO
rule chrom_info_detail:
    input:
        gtf = config['ref']['annotation']
    output:
        gtf = annotation_dir + "/chrom_info_detail.xls"
    run:
        cmd = "python scripts/diff/construct_chrom_info.py {input.gtf} {output.gtf}"
        shell(cmd)


rule ggbio_circle:
    input:
        fpkm = exp_dir + "/gene_fpkm_all_samples_anno.txt",
        metadata = config["metadata"],
        chrom_info = config["ref"]["chrom_info"],
        chrom_info_detail = annotation_dir + "/chrom_info_detail.xls",
        deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_anno.xls",
    params:
        contrast = "{contrast}",
        ggbio_parameter = config['ggbio_parameter'],
    output:
        circle_pdf = diff_dir + "/{level}/{method}/{contrast}/{contrast}_circle.pdf",
        circle_png = diff_dir + "/{level}/{method}/{contrast}/{contrast}_circle.png",
    shell:
        "/nfs1/public2/User/axl/miniconda3/bin/Rscript scripts/diff/ggbio.R {input.fpkm} {input.metadata} {input.chrom_info} {input.chrom_info_detail} {params.contrast} {input.deg} {params.ggbio_parameter} {output.circle_pdf} {output.circle_png} "


rule cog_plot:
    input:
        cog_anno = annotation_dir + "/COG/COG_annotations.xls",
        deg_list = diff_dir + "/{level}/{method}/DEGs_list/all_DEGs_list.xls",
    output:
        cog_dir = directory(diff_dir + '/{level}/{method}/{contrast}/cog/'),
    params:
        contrast = "{contrast}",
    shell:
        "/nfs1/public2/User/wzc/miniconda3/envs/rnaseq/bin/Rscript scripts/diff/cog_anno_barplot.R -c {params.contrast} -a {input.cog_anno} -d {input.deg_list} -o {output.cog_dir}"



# new_add
# rule PPI:
#     input:
#         deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
#         fasta = annotation_dir + "/genes.fa",
#     params:
#         species_code = config['species_code'],
#         contrast = '{contrast}'
#     output:
#         temp = temp(diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_gene_lst.xls"),
#         PPI_res = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_PPI.xls"
#     run:
#         cwd = os.getcwd()
#         cmd = """awk -F '\\t' '{{if($9 == "Decreased" || $9 == "Increased") print $8}}' {input.deg} > {output.temp} && """
#         cmd += "/nfs1/public2/Pipe2/miniconda3/envs/noref-ppi/bin/python scripts/diff/BLASTX_TO_PPI_v5.11.py --entries %s/{output.temp} --species {params.species_code} --fa %s/{input.fasta} --name %s/{params.contrast} --output %s/{output.PPI_res}" % (cwd, cwd,cwd, cwd)
#         shell(cmd)

# rule TF:
#     input:
#         deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_all.xls",
#         fasta = annotation_dir + "/genes.fa",
#         target = config["ref"]["tftar"],
#         TF_list = config["ref"]["TF_list"],
#     params:
#         contrast = '{contrast}'
#     output:
#         temp = temp(diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_gene_lst.xls"),
#         TF_res = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_TF.xls"
#     run:
#         cwd = os.getcwd()
#         cmd = """awk -F '\\t' '{{if($9 == "Decreased" || $9 == "Increased") print $8}}' {input.deg} > {output.temp} && """
#         cmd += "/nfs1/public2/Pipe2/miniconda3/envs/noref-ppi/bin/python scripts/diff/BLASTX_TO_TF_v5.11.py --entries %s/{output.temp} --fa %s/{input.fasta} --TF {input.TF_list} --name %s/{params.contrast} --target {input.target} --output %s/{output.TF_res}" % (cwd, cwd,cwd, cwd)
#         shell(cmd)



# rule add_TF_annotation_diff:
#     input:
#         deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
#         anno = config["ref"]["TF_list"]
#     output:
#         diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_TF.xls",
#     params:
#         Rscript = config["software"]["Rscript"]
#     shell:
#         "Rscript scripts/diff/add_tf_annotation.R {input.deg} {input.anno} {output}"

# rule PPI:
#     input:
#         deg = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff.xls",
#         PPI_anno = config['ref']['PPIanno']
#     params:
#         score_limit = config['PPI_score'], 
#         label = '{label}',
#         python = config["software"]["python"]
#     output:
#         PPI_table = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_PPI_{label}.xls"
#     shell:
#         "{params.python} scripts/diff/merge_PPI.py {input.deg} {input.PPI_anno} {params.label} {params.score_limit} {output.PPI_table}"


