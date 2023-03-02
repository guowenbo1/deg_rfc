
#kegg通路富集图
import time

# ruleorder:kegg_enrich > weblink_anno


rule kegg_top_20:
    input:
        diff_dir +"/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_all.xls"
    output:
        temp(diff_dir +"/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_all_draw_no_submit.xls")
    shell:
        "head -21 {input} > {output}"


rule weblink_anno:
    input:
        gene2pathway = annotation_dir + "/kegg/kegg_anno_gene2pathway.txt",
        kegg_file = diff_dir + "/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_all.xls",
        diff_anno = diff_dir + "/{level}/{method}/{contrast}/{contrast}_diff_anno.xls",
    params:
        python = config["software"]["python"]
    output:
        output1 = diff_dir + "/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_all_weblink.xls",
    shell:    
        "python scripts/diff/autoAddweblinkOfkegg.py -a {input.gene2pathway} -k {input.kegg_file} -d {input.diff_anno}"        



rule kegg_mapplot:
    input:
        diff_anno = diff_dir +"/{level}/{method}/{contrast}/{contrast}_diff_anno.xls",
        kegg_file = diff_dir +"/{level}/{method}/{contrast}/kegg_enrich/kegg_enrich_all_draw_no_submit.xls",
    output:
        directory(diff_dir +"/{level}/{method}/{contrast}/kegg_enrich/maps/")          
    params:
        python = config["software"]["python"]
    run:        
        cmd = "python scripts/diff/KEGGAnnoAuto.py {input.diff_anno} {input.kegg_file} {output}"  
        shell(cmd)  
