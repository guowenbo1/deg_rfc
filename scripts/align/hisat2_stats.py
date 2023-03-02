#!/usr/bin/env python
import os, sys, re

# inputdir, out = sys.argv[1:]
files = sys.argv[1:len(sys.argv)]
out = sys.argv[-1]
# inputdir = "."
# out = "out.txt"

def get_hisat2_stats(file):
    with open(file, "r") as f:
        data = []
        for line in f:
            line = line.strip()
            num = line.split()[0]
            data.append(num)
    total = int(data[0])*2
    total_ratio = data[-1]
    uniq_ratio = (int(data[3])*2 + int(data[7])*2 + int(data[12]))/total
    sample = os.path.basename(file).split(".")[0]
    return {"sample":sample,"total":total,"total_ratio":total_ratio,"uniq_ratio":uniq_ratio}

# files = os.listdir(inputdir)
with open(out, "w") as f:
    print("Sample\tTotal Clean Reads\tTotal Mapping Ratio\tUniquely Mapping Ratio", file = f)
    for file in sorted(files):
        if file.endswith(".log"):
            stats_dict = get_hisat2_stats(file)
            f.write(f"{stats_dict['sample']}\t{stats_dict['total']}\t{stats_dict['total_ratio']}\t{stats_dict['uniq_ratio']*100:.2f}%\n")

