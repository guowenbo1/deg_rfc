# new_trans

rule get_new_trans:
    input:
        gtf_dir = assembly_dir + '/comparison/'
    output:
        new_trans = new_trans_dir + "/new_trans.gtf",
    params:
        tmap = lambda wildcards, input: assembly_dir + "/merged/all.merged.gtf.tmap",
    shell:
        "python {workflow.basedir}/scripts/new_trans/filter_new_trans.py  {params.tmap} {input.gtf_dir}/all.annotated.gtf {output} "