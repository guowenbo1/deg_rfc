# map trimmed fastq to genome

rule align:
    input:
        # fastq1 = clean_data_dir + "/{sample}_R1.fq.gz",
        # fastq2 = clean_data_dir + "/{sample}_R2.fq.gz",
        fastq1 = clean_data_dir + "/fastp/{sample}_R1.fq.gz",
        fastq2 = clean_data_dir + "/fastp/{sample}_R2.fq.gz",
    output:
        align_dir + "/{sample}.bam",
    log:
        align_dir + "/hisat2/{sample}.log",
    params:
        library_type = config["library_type"],
        # ref = config["ref"]["index"],
        threads = config['threads'],
        hisat2 = config["software"]["hisat2"],
        samtools = config["software"]["samtools"],
    shell:
        '''{params.hisat2} -p {params.threads} --dta --rna-strandness {params.library_type} -x {params.ref} \
        -1 {input.fastq1} -2 {input.fastq2} --summary-file {log} | {params.samtools} sort -@ {params.threads} -o {output}'''


rule align_stats:
    input:
        expand(align_dir + "/hisat2/{sample}.log", sample=samples),
    output:
        align_dir + "/hisat2/hisat2_stats.xls",
    params:
        python = config["software"]["python"],
    shell:
        "{params.python} {workflow.basedir}/scripts/align/hisat2_stats.py {input} {output} "


rule gtf2bed:
    input:
        gtf = config["ref"]["annotation"],
    output:
        bed = align_dir + "/QC/ref.bed",
    params:
        perl = config["software"]["perl"],
    shell:
        "{params.perl} {workflow.basedir}/scripts/align/gtf2bed.pl {input.gtf} > {output.bed}"


# 未测试
rule samtools_index:
    input:
        align_dir + "/{sample}.bam",
    output:
        align_dir + "/{sample}.bam.bai",
    params:
        threads = config["threads"],
        samtools = config["software"]["samtools"]
    run:
        shell("{params.samtools} index -@ {params.threads} {input}")
        csi = "{input}.csi"
        cmd = "mv" + csi + "{output}"
        if os.path.exists(csi):
            shell(cmd)
# rule samtools_index:
#     input:
#         align_dir + "/{sample}.bam"
#     output:
#         csi = align_dir + "/{sample}.bam.csi",
#         bai = align_dir + "/{sample}.bam.bai",
#     params:
#         threads = config["threads"]
#     shell:
#         "samtools index -@ {params.threads} {input} -c {output.csi} && mv {output.csi} {output.bai}"

rule junctionSaturation:
    input:
        bam = align_dir + "/{sample}.bam",
        bed = align_dir + "/QC/ref.bed"
    output:
        rlan = align_dir + "/QC/junctionSaturation/{sample}.junctionSaturation_plot.r",
    run:
        prefix = output.rlan.replace(".junctionSaturation_plot.r", "")
        cmd = "junction_saturation.py -i {input.bam} -r {input.bed} -o {prefix}"
        shell(cmd)
        # pdf2png
        # cmd = "pdftoppm {output.res} {prefix} -png -f 1 && mv {prefix}-1.png {output.png}"
        # shell(cmd)

rule junction_pic:
    input:
        rlan = align_dir + "/QC/junctionSaturation/{sample}.junctionSaturation_plot.r",
    output:
        res = align_dir + "/QC/junctionSaturation/{sample}.junctionSaturation_plot.pdf",
        png = align_dir + "/QC/junctionSaturation/{sample}.junctionSaturation_plot.png",
    run:
        prefix = output.res.replace(".junctionSaturation_plot.pdf", "")
        cmd = "Rscript {input.rlan}"
        shell(cmd)
        cmd = "/nfs1/public2/User/JXQ/miniconda3/bin/pdftoppm {output.res} {prefix} -png -f 1 && mv {prefix}-1.png {output.png}"
        shell(cmd)

rule geneBody_coverage:
    input:
        bam = expand(align_dir + "/{sample}.bam", sample = samples),
        bai = expand(align_dir + "/{sample}.bam.bai", sample = samples),
        bed = align_dir + "/QC/ref.bed"
    output:
        expand(expand(align_dir + "/QC/geneBody_coverage/{sample}.geneBodyCoverage.curves.{{exts}}", sample=samples), exts = fig_exts)
    shell:
        "bash ./scripts/align/gene_body.sh"

rule read_duplication:
    input:
        bam = align_dir + "/{sample}.bam",
        # bed = align_dir + "QC/ref.bed",
    output:
        # pos = align_dir + "/QC/read_duplication/{sample}.pos.DupRate.xls",
        # seq = align_dir + "/QC/read_duplication/{sample}.seq.DupRate.xls",
        rlan = align_dir + "/QC/read_duplication/{sample}.DupRate_plot.r",
    run:
        prefix = output.rlan.replace(".DupRate_plot.r", "")
        cmd  = "read_duplication.py -i {input.bam} -o {prefix}"
        shell(cmd)
        # pdf2png
        # cmd = "pdftoppm {output.pdf} {prefix} -png -f 1 && mv {prefix}-1.png {output.png}"
        # shell(cmd)

rule dupli_pic:
    input:
        rlan = align_dir + "/QC/read_duplication/{sample}.DupRate_plot.r",
    output:
        pdf = align_dir + "/QC/read_duplication/{sample}.DupRate_plot.pdf",
        png = align_dir + "/QC/read_duplication/{sample}.DupRate_plot.png",
    run:
        prefix = output.pdf.replace(".DupRate_plot.pdf", "")
        cmd = "Rscript {input.rlan}"
        shell(cmd)
        cmd = "/nfs1/public2/User/JXQ/miniconda3/bin/pdftoppm {output.pdf} {prefix} -png -f 1 && mv {prefix}-1.png {output.png}"
        shell(cmd)


rule RPKM_saturation:
    input:
        bam = align_dir + "/{sample}.bam",
        bed = align_dir + "/QC/ref.bed"
    output:
        eRPKM = align_dir + "/QC/RPKM_saturation/{sample}.eRPKM.xls",
        rawCount = align_dir + "/QC/RPKM_saturation/{sample}.rawCount.xls",
        rlan = align_dir + "/QC/RPKM_saturation/{sample}.saturation.r",
    run:
        prefix = output.rlan.replace(".saturation.r", "")
        cmd = "RPKM_saturation.py -r {input.bed} -i {input.bam} -o {prefix}"
        shell(cmd)
        # pdf2png
        # cmd = "pdftoppm {output.pdf} {prefix} -png -f 1 && mv {prefix}-1.png {output.png}"
        # shell(cmd)

rule RPKM_pic:
    input:
        rlan = align_dir + "/QC/RPKM_saturation/{sample}.saturation.r",
    output:
        pdf = align_dir + "/QC/RPKM_saturation/{sample}.saturation.pdf",
        png = align_dir + "/QC/RPKM_saturation/{sample}.saturation.png",
    run:
        prefix = output.pdf.replace(".saturation.pdf", "")
        cmd = "Rscript {input.rlan}"
        shell(cmd)
        cmd = "/nfs1/public2/User/JXQ/miniconda3/bin/pdftoppm {output.pdf} {prefix} -png -f 1 && mv {prefix}-1.png {output.png}"
        shell(cmd)


rule RPKM_saturation_plot:
    input:
        rawCount = align_dir + "/QC/RPKM_saturation/{sample}.rawCount.xls",
    output:
        fig = expand(align_dir + "/QC/RPKM_saturation/{{sample}}_saturation_curve.{ext}", ext = fig_exts)
    script:
        "../scripts/align/RPKM_saturation_plot.R"

rule read_distribution:
    input:
        bam = align_dir + "/{sample}.bam",
        bed = align_dir + "/QC/ref.bed",
    output:
        align_dir + "/QC/read_distribution/{sample}.read_distribution.txt"
    shell:
        "read_distribution.py  -i {input.bam} -r {input.bed} > {output}"

rule read_distribution_stat:
    input:
        align_dir + "/QC/read_distribution/{sample}.read_distribution.txt"
    output:
        align_dir + "/QC/read_distribution/{sample}.read_distribution_stat.txt"
    shell:
        "perl scripts/align/read_distribution_stat.pl  {input} {output}"

rule read_distribution_plot:
    input:
        align_dir + "/QC/read_distribution/{sample}.read_distribution_stat.txt"
    output:
        pdf = align_dir + "/QC/read_distribution/{sample}.read_distribution_plot.pdf",
        png = align_dir + "/QC/read_distribution/{sample}.read_distribution_plot.png"
    shell:
        "Rscript scripts/align/read_distribution.R  {input} {output.pdf} {output.png} "
