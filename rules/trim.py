def get_fastq(wildcards):
    return metadata.loc[wildcards.sample, ["fq1", "fq2"]].dropna()


rule fastq_stats:
    input:
        get_fastq
    output:
        seqkit_stats_dir + "/{sample}.xls"
    params:
        threads = 10
    shell:
        " cat {input}| seqkit stats -aTj {params.threads} -| sed 's/-/{wildcards.sample}/' -> {output} "


# rule data_enough:
#     input:
#         seqkit_stats_dir + "/{sample}.xls"
#     output:

rule fastp_pe:
    input:
        get_fastq
    output:
        fastq1 = clean_data_dir + "/fastp/{sample}_R1.fq.gz",
        fastq2 = clean_data_dir + "/fastp/{sample}_R2.fq.gz",
        html = clean_data_dir + "/fastp/html/{sample}.html",
        json = clean_data_dir + "/fastp/json/{sample}.json",
    params:
        adapter = config["adapter"],
        truc_f = config["truc_f"],
        fastp = config["software"]["fastp"],
        truc_F = 0,
        len_min = 50,
        threads = config['threads'],
    shell:
        """{params.fastp} -i {input[0]} -I {input[1]} -a {params.adapter} \
            -f {params.truc_f} -F {params.truc_F} -n 5 -l {params.len_min} \
            -3 -W 4 -M 20 -q 20 -w {params.threads} -h {output.html} -j {output.json} \
            -o {output.fastq1} -O {output.fastq2}"""


rule derRNA:
    input:
        fastq1 = clean_data_dir + "/fastp/{sample}_R1.fq.gz",
        fastq2 = clean_data_dir + "/fastp/{sample}_R2.fq.gz",
    output:
        fastq1 = clean_data_dir + "/{sample}_R1.fq.gz",
        fastq2 = clean_data_dir + "/{sample}_R2.fq.gz",
    params:
        rrnadb = config["rrnadb"],
        threads = config["threads"]
    shell:
        ''' bowtie2 -p {params.threads} -x {params.rrnadb} -1 {input.fastq1} -2 {input.fastq2} \
            -S /dev/null --un-conc-gz {wildcards.sample} && \
            mv {wildcards.sample}.1 {output.fastq1} && mv {wildcards.sample}.2 {output.fastq2}'''


rule fastp_stats:
    input:
        json = expand(clean_data_dir + "/fastp/json/{sample}.json", sample=samples)
    output:
        clean_data_dir + "/fastp/fastp_stats.xls"
    shell:
        " python {workflow.basedir}/scripts/trim/fastp_stats.py {input.json} {output} "


rule derRNA_seqkit_stats:
    input:
        fastq1 = clean_data_dir + "/{sample}_R1.fq.gz",
        fastq2 = clean_data_dir + "/{sample}_R2.fq.gz",
    output:
        clean_data_dir + "/stats/{sample}.xls"
    params:
        threads = config["threads"]
    shell:
        " cat {input}| seqkit stats -aTj {params.threads} -| sed 's/-/{wildcards.sample}/' -> {output} "


rule derRNA_stats_summary:
    input:
        expand(clean_data_dir + "/stats/{sample}.xls", sample=samples)
    output:
        clean_data_dir + "/derRNA_stats.xls"
    run:
        with open(output[0], "w") as f:
            line = "sample\tderRNA_data_reads\tderRNA_data_bases\tderRNA_data_q20_rate\tderRNA_data_q30_rate\n"
            f.write(line)
            for file in input:
                with open(file) as f_obj:
                    for line in f_obj:
                        if line.startswith("file"):
                            continue
                        line_lst = line.split("\t")
                        line = "\t".join([line_lst[0], line_lst[3], "{:.2f} G".format(int(line_lst[4]) / 1_000_000_000), line_lst[13], line_lst[14]])
                        f.write(line)


rule filter_stats_summary:
    input:
        clean_data_dir + "/fastp/fastp_stats.xls",
        clean_data_dir + "/derRNA_stats.xls"
    output:
        clean_data_dir + "/filter_stats.xls"
    run:
        df1 = pd.read_csv(input[0], sep="\t")
        df2 = pd.read_csv(input[1], sep="\t")
        df = pd.merge(df1, df2, on="sample")
        df.to_csv(output[0], index=False, sep="\t")

rule qc_stat:
    input:
        json = clean_data_dir + "/fastp/json/{sample}.json"
    params:
        sample = "{sample}"
    output:
        c_GC = temp(clean_data_dir + "/fastp/json/clean_{sample}.GC"),
        c_QM = temp(clean_data_dir + "/fastp/json/clean_{sample}.QM"),
        r_GC = temp(clean_data_dir + "/fastp/json/raw_{sample}.GC"),
        r_QM = temp(clean_data_dir + "/fastp/json/raw_{sample}.QM"),
        stat = temp(clean_data_dir + "/fastp/json/{sample}.stat"),
    shell:
        "/nfs1/public2/Pipe2/miniconda3/envs/noref-qc/bin/python /nfs2/public2/User/lixiaofei/projects/new-mRNA/mRNA/bin/fastp_convert --json {input.json} --identity {params.sample}"

rule qc_stat_pic:
    input:
        c_GC = clean_data_dir + "/fastp/json/clean_{sample}.GC",
        c_QM = clean_data_dir + "/fastp/json/clean_{sample}.QM",
        r_GC = clean_data_dir + "/fastp/json/raw_{sample}.GC",
        r_QM = clean_data_dir + "/fastp/json/raw_{sample}.QM",
        stat = clean_data_dir + "/fastp/json/{sample}.stat",
    params:
        sample = "{sample}",
        outdir = clean_data_dir + "/fastp/{sample}/"
    output:
        Error_png = clean_data_dir + "/fastp/{sample}/{sample}.Error.png",
        GC_png = clean_data_dir + "/fastp/{sample}/{sample}.GC.png",
    shell:
        "/nfs2/public2/User/lixiaofei/projects/new-mRNA/mRNA/bin/ngqc_plot --gc {input.c_GC} --qm {input.c_QM} --stat {input.stat} --name {params.sample} --outdir {params.outdir}"



