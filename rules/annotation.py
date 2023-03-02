# function annotation

# NR NT  uniprot COG KEGG GO

databases = ["NR", "Swiss-Prot", "TrEMBL", "EggNOG"]


rule get_trans_genes:
    input:
        gtf = config['ref']['annotation'],
        genome = config['ref']['genome']
    output:
        trans = annotation_dir + "/transcripts.fa",
        gene = annotation_dir + "/genes.fa"
    shell:
        "gffread -g {input.genome} -w {output.trans} {input.gtf} &&"
        "python {workflow.basedir}/scripts/annotation/get_genes.py {output.trans} {output.gene}"


rule get_gene_id_name:
    input:
        gtf = config["ref"]["annotation"],
    output:
        annotation_dir + "/gene_id2name.txt"
    shell:
        "python {workflow.basedir}/scripts/annotation/get_id_new.py {input.gtf} {output}"


rule eggnog:
    input:
        gene = annotation_dir + "/genes.fa"
    output:
        annotation_dir + "/eggnog/eggnog.emapper.annotations",
        annotation_dir + "/eggnog/eggnog.emapper.seed_orthologs",
    params:
        prefix = "eggnog",
        threads = config["threads"]
    run:
        outdir = os.path.dirname(output[0])
        shell("emapper.py --translate -i {input.gene} --cpu {params.threads} -m diamond -o {params.prefix} --output_dir {outdir} ")
        # shell("emapper.py -i {input.gene} --cpu {params.threads} -m diamond -o {params.prefix} --output_dir {outdir} ")


rule eggnog_format:
    input:
        annotation_dir + "/eggnog/eggnog.emapper.annotations",
        annotation_dir + "/gene_id2name.txt"
    output:
        annotation_dir + "/eggnog/eggnog.txt",
    run:
        f1=open(input[0])
        f2=open(output[0], "w")
        line = "gene_name,seed_eggNOG_ortholog,seed_ortholog_evalue,seed_ortholog_score,best_tax_level,"
        line+="Preferred_name,GOs,EC,KEGG_ko,KEGG_Pathway,KEGG_Module,KEGG_Reaction,KEGG_rclass,BRITE,KEGG_TC,"
        line+="CAZy,BiGG_Reaction,tax_scope,OGs,bestOG,COG_cat,eggNOG_annot"
        line = "\t".join(line.split(",")) + "\n"
        f2.write(line)
        for line in f1:
            if line.startswith("#"):
                continue
            f2.write(line)
        f1.close()
        f2.close()
        df1 = pd.read_csv(output[0], sep='\t')
        df2 = pd.read_csv(input[1], sep='\t')
        df = pd.merge(df2, df1, how="left", on="gene_name")
        df.to_csv(output[0], index=False, sep='\t')


rule go_anno:
    input:
        eggnog = annotation_dir + "/eggnog/eggnog.txt",
    output:
        go_anno = annotation_dir + "/go/go.txt",
        go_anno_level2 = annotation_dir + "/go/go_level2.txt",
        go2gene = annotation_dir + "/go/go2gene_level2.txt",
        fig = expand(annotation_dir + "/go/go_level2_barplot.{ext}", ext=fig_exts)
    script:
        "../scripts/annotation/go_anno.R"


rule gene2go:
    input:
        go_anno = annotation_dir + "/go/go.txt",
    output:
        gene2go = annotation_dir + "/go/gene2go.xls",
    shell:
        '''csvtk mutate2 -t -l -j 20 -n GO -e '$GOID + "," + $ONTOLOGY + "," + $TERM' {input.go_anno}'''
        '''|csvtk cut -t -j 20 -f gene_id,gene_name,GO|csvtk collapse -t -j 20 -f gene_id,gene_name -v GO > {output.gene2go}'''

rule kegg_anno:
    input:
        eggnog = annotation_dir + "/eggnog/eggnog.txt",
        kegg = config["kegg_db"]
    output:
        kegg_anno = annotation_dir + "/kegg/kegg.txt",
        gene2pathway = annotation_dir + "/kegg/kegg_anno_gene2pathway.txt",
        type2gene = annotation_dir + "/kegg/kegg_anno_type2gene.txt",
        pathway2gene = annotation_dir + "/kegg/kegg_anno_pathway2gen2.txt",
        fig = expand(annotation_dir + "/kegg/kegg_level2_barplot.{ext}", ext=fig_exts),
    script:
        "../scripts/annotation/kegg_anno.R"


rule gene2kegg:
    input:
        gene2pathway = annotation_dir + "/kegg/kegg_anno_gene2pathway.txt",
    output:
        gene2kegg = annotation_dir + "/kegg/gene2kegg.xls",
    shell:
        '''csvtk mutate2 -t -l -j 20 -n KO -e '$ko + "," + $ko_name + "," + $ko_des' {input}|csvtk collapse -t -j 20 -f gene_id,gene_name -v KO,pathway > {output}'''


rule get_all_annotation:
    input:
        gene2go = annotation_dir + "/go/gene2go.xls",
        gene2kegg = annotation_dir + "/kegg/gene2kegg.xls",
    output:
        anno = annotation_dir + "/gene2function.xls",
        anno_temp = temp(annotation_dir + "/gene2function_temp.xls")
    shell:
        "/nfs1/public2/User/wzc/apps/bin/csvtk join -t -l {input} > {output.anno_temp} && Rscript ./scripts/annotation/get_all_annotation.R {input} {output.anno}"


rule cog_anno_all:
    input:
        gene = annotation_dir + "/genes.fa",
    output:
        cog_anno = annotation_dir + "/COG/COG_annotations.xls",
    shell:
        "python {workflow.basedir}/scripts/annotation/run_cog.py -q {input.gene} -t fna -vp Annotation/COG/COG_anno_all -o {output.cog_anno}"
