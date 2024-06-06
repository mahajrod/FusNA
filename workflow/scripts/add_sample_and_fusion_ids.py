#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys, gzip, bz2

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file
import argparse
import pandas as pd


def metaopen(filename, flags, buffering=None, compresslevel=5):
    if not isinstance(filename, str): # or isinstance(filename, gzip.GzipFile) or isinstance(filename, bz2.BZ2File):
        if isinstance(filename, file):
            return filename
        else:
            raise ValueError("ERROR!!! Not file object or str: {}".format(str(filename)))
    elif filename[-3:] == ".gz":
        return gzip.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    elif filename[-4:] == ".bz2":
        return bz2.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    else:
        if buffering is not None:
            return open(filename, flags, buffering=buffering)
        else:
            return open(filename, flags)


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="File with fusions called by Arriba. Default: stdin")
parser.add_argument("-s", "--sample_id", action="store", dest="sample_id", required=True,
                    help="Sample id to add to the fusions as a first column")
parser.add_argument("-r", "--fusion_id_suffix", action="store", dest="fusion_id_suffix", default=".FUS",
                    help="Fusion id suffix (prefix is sample id). Default: '.FUS' ")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output file. Default: stdout")
args = parser.parse_args()

fusion_index = 1
with metaopen(args.input, "r") as in_fd, metaopen(args.output, "w") as out_fd:
    header_line = in_fd.readline()
    out_fd.write(f"#sample_id\tfusion_id\t{header_line[1:]}")
    for line in in_fd:
        out_fd.write(f"{args.sample_id}\t{args.sample_id}{args.fusion_id_suffix}{fusion_index}\t" + line)
        fusion_index += 1
