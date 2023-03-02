import sys



tmap, gtf, out = sys.argv[1:]

classcode = ["u", "x", "i", "j", "o"] # ixuoe

list_qry_id = []

with open(tmap) as file_obj:

    i = 0

    for line in file_obj:

        i += 1

        if i >= 2:

            splits = line.split("\t")

            # if splits[2] in classcode and int(splits[5]) >= 2 and float(splits[6]) >= 0.5 and int(splits[9]) > 200:

            if splits[2] in classcode and int(splits[5]) >= 2 and int(splits[9]) > 200:

                list_qry_id.append(splits[4])



with open(gtf) as file_obj:

    with open(out, "w") as file_rs:

        for line in file_obj:

            if not line.startswith("#"):

                splits = line.split("\t")

                trans_id_seq = splits[8].split(";")

                transcript_id = trans_id_seq[0].split("\"")[1]

                if transcript_id in list_qry_id:

                    file_rs.write(line)
