#!/usr/bin/env python

"""
Read a directory of .hla files, clean them by stripping them of
extraneous information (apostrophes and numerical values), and output
the updated .hla files in a new directory.

<this script> --in <input_dir> --out <output_dir>
"""

import argparse
from os import listdir, path

def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--indir', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    for filename in listdir(args.indir):
        file_path = path.join(args.indir, filename)
        lines = [strip_info(line) for line in open(file_path)]
        new_path = path.join(args.outdir, filename)
        with open(new_path, 'a') as nf:
            for line in lines:
                nf.write(line)

def strip_info(line):
    line = line.strip()
    line = line.split(' ')[0]
    line = line.split('\'')[0]
    return line + '\n'

if __name__ == "__main__":
    run()

