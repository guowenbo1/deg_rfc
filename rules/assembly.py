# assembly

rule assembly:
    input:
        align_dir + "/{sample}.bam",
    output:
        assembly_dir + "/{sample}/transcripts.gtf",
    params:
        gtf = config["ref"]["annotation"],
        string_type = config["string_type"],
        threads = config['threads'],
        stringtie = config["software"]["stringtie"],
    shell:
        "{params.stringtie} -p {params.threads} {params.string_type} {params.gtf} -l {wildcards.sample} -o {output} {input}"


rule compose_merge:
    input:
        expand(assembly_dir + "/{sample}/transcripts.gtf", sample=samples),
    output:
        txt = assembly_dir + '/assemblies.txt',
    run:
        with open(output.txt, 'w') as out:
            # print到文件
            print(*input, sep="\n", file=out)


# stringtie --merge defaults -m <min_len> 50;-c <min_cov> 0;-F <min_fpkm> 1;-T <min_tpm> 1;
rule merge_assemblies:
    input:
        assembly_dir + '/assemblies.txt',
    output:
        assembly_dir + '/merged/merged.gtf',
    params:
        length = 200,
        gtf = config["ref"]["annotation"],
        threads = config['threads'],
        stringtie = config["software"]["stringtie"],
    shell:
        "{params.stringtie} --merge -p {params.threads} -m {params.length} -G {params.gtf} -o {output} {input}"


# gffcomapre replace cuffcompare
rule compare_assemblies:
    input:
        assembly_dir + '/merged/merged.gtf',
    output:
        directory(assembly_dir + '/comparison/')
    params:
        gtf = config["ref"]["annotation"],
        genome = config["ref"]["genome"],
        gffcompare = config["software"]["gffcompare"],
    shell:
        'mkdir -p {output} && {params.gffcompare} -o {output}/all -s {params.genome} -r {params.gtf} {input}'


rule stringtie_expression:
    input:
        config["ref"]["annotation"],
        align_dir + "/{sample}.bam",
    output:
        assembly_dir + "/ballgrown/{sample}/transcripts.gtf",
        assembly_dir + "/ballgrown/{sample}/gene_abundances.tsv",
    params:
        threads = config["threads"],
        stringtie = config["software"]["stringtie"]
    shell:
        "{params.stringtie} -e -B -p {params.threads} -G {input[0]} -A {output[1]} -o {output[0]} {input[1]} "
