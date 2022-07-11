import getopt
import pandas as pd
import sys
from io import StringIO

class myGffParser():
    def __init__(self):
        input_file = sys.argv[1:]
        opts, _ = getopt.getopt(input_file, "i:", ['input'])
        print(opts)
        for opt, arg in opts:
            if opt in ["-i", "--input"]:
                with open(arg) as f:
                    lines = f.readlines()
        df = pd.DataFrame(columns=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

        for i, line in enumerate(lines):
            if line[0] != "#" and len(line.strip().split("\t")) == 8:
                record = line.strip().split("\t")
                df.loc[i, "seqid"] = record[0]
                df.loc[i, "source"] = record[1]
                df.loc[i, "type"] = record[2]
                df.loc[i, "start"] = int(record[3])
                df.loc[i, "end"] = int(record[4])
                df.loc[i, "score"] = record[5]
                df.loc[i, "strand"] = record[6]
                df.loc[i, "phase"] = record[7][0]
                df.loc[i, "attributes"] = record[7][1:]
        self.df = df
        print(df)