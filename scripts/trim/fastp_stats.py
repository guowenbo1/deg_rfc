#!/usr/bin/env python
import sys
import json
import os

files = sys.argv[1:len(sys.argv)-1]
outfile = sys.argv[-1]

head = "sample\traw_data_reads\traw_data_bases\traw_data_q20_bases\traw_data_q30_bases\traw_data_q20_rate\traw_data_q30_rate"
head += "\tclean_data_reads\tclean_data_bases\tclean_data_q20_bases\tclean_data_q30_bases\tclean_data_q20_rate\tclean_data_q30_rate\n"

with open(outfile, "w") as of:
    of.write(head)
    for file in files:
        # print(file)
        sample = os.path.basename(file).replace(".json", "")
        with open(file) as jh:
            data = json.load(jh)
        b1 = data["summary"]["before_filtering"]["total_reads"]
        b2 = data["summary"]["before_filtering"]["total_bases"]
        b2 = "{:.2f}".format(b2 / 1_000_000_000) + "G"
        b3 = data["summary"]["before_filtering"]["q20_bases"]
        b4 = data["summary"]["before_filtering"]["q30_bases"]
        b5 = data["summary"]["before_filtering"]["q20_rate"]
        b6 = data["summary"]["before_filtering"]["q30_rate"]
        a1 = data["summary"]["after_filtering"]["total_reads"]
        a2 = data["summary"]["after_filtering"]["total_bases"]
        a2 = "{:.2f}".format(a2 / 1_000_000_000) + "G"
        a3 = data["summary"]["after_filtering"]["q20_bases"]
        a4 = data["summary"]["after_filtering"]["q30_bases"]
        a5 = data["summary"]["after_filtering"]["q20_rate"]
        a6 = data["summary"]["after_filtering"]["q30_rate"]
        s = sample +"\t" + "\t".join(map(str, [b1, b2, b3, b4, b5, b6, a1, a2, a3, a4, a5, a6]))
        s += "\n"
        of.write(s)
