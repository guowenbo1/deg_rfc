rule get_bam_lst:
    input:
        metadata = config["metadata"],
        diff_lst = config["diff_lst"]
    output:
        b1 = as_dir + "/{contrast}/{contrast}_b1.txt",
        b2 = as_dir + "/{contrast}/{contrast}_b2.txt"
    run:
        diff_lst = input.diff_lst
        diff_groups = []
        contrasts = []
        contrasts_groups = {}
        with open(diff_lst, "r") as f:
            for line in f:
                line_list = line.split()
                contrasts.append(line_list[0])
                diff_groups.append(line_list[1])
                contrasts_groups[line_list[0]] = [line_list[2], line_list[3]]
        group_type = diff_groups[contrasts.index(wildcards.contrast)]
        groups = contrasts_groups[wildcards.contrast]
        group1 = groups[0]
        group2 = groups[1]
        metadata = pd.read_csv(input.metadata, sep="\t").set_index("SampleID", drop=False)
        metadata1 = metadata[metadata[group_type]==group1]
        metadata2 = metadata[metadata[group_type]==group2]
        f1 = open(output.b1, "w")
        f2 = open(output.b2, "w")
        samples1 = []
        samples2 = []
        for sample in metadata1.index:
            line = align_dir + "/" + sample + ".bam"
            samples1.append(line)
        f1.write(",".join(samples1))
        f1.close()
        for sample in metadata2.index:
            line = align_dir + "/" + sample + ".bam" 
            samples2.append(line)
        f2.write(",".join(samples2))
        f2.close()


rule rmats:
    input:
        bam = expand(align_dir + "/{sample}.bam", sample=samples),
        b1 = as_dir + "/{contrast}/{contrast}_b1.txt",
        b2 = as_dir + "/{contrast}/{contrast}_b2.txt",
        gtf = config['ref']['annotation'],
    output:
        directory(as_dir + "/{contrast}/results")
    params:
        libType = config["libType"],
        threads = config["threads"],
        tmp = directory(as_dir + "/{contrast}/tmp")
    shell:
         # "export PATH=/nfs1/public2/User/JXQ/mRNA-seq/HT20200717084-Yangmiqi/.snakemake/conda/37717c69:$PATH &&"
         "/nfs1/public2/User/felix/biosoft/miniconda3/envs/rmats/bin/python /nfs1/public2/User/felix/biosoft/miniconda3/envs/rmats/bin/rmats.py --b1 {input.b1} --b2 {input.b2} --gtf {input.gtf} --od {output} -t paired --novelSS --nthread {params.threads} --readLength 142 --variable-read-length --tmp {params.tmp}"

         # "/nfs1/public2/User/JXQ/mRNA-seq/HT20200717084-Yangmiqi/.snakemake/conda/37717c69/bin/rmats.py "
