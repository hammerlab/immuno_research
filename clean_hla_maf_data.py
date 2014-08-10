#!/usr/bin/env python

"""
Reads HLA and MAF file directories, and does the following:

    * Cleans HLA files by stripping them of extraneous information
    * Removes any HLA files without a corresponding MAF, and vice versa
    * Converts .maf.txt to .maf
    * Strip HLaA filenames down to 15 chars + extension
    * Writes new HLA and MAF files to new directories with the suffix
      "_cleaned"

<this script> --hla-dir <hla-dir> --maf-dir <maf-dir>
"""

import argparse
from os import listdir, path, rename, makedirs, extsep
import shutil

# The number of characters (excluding extension) to keep for HLA files
NUM_CHARS_HLA = 15

# TODO Make this code less ugly. Separate into different functions, etc.
def run():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument('--hla-dir', required=True)
    parser.add_argument('--maf-dir', required=True)
    args = parser.parse_args()
    hla_dir_out = args.hla_dir + "_cleaned"
    maf_dir_out = args.maf_dir + "_cleaned"
    shutil.rmtree(hla_dir_out, ignore_errors=True)
    shutil.rmtree(maf_dir_out, ignore_errors=True)
    makedirs(hla_dir_out)
    makedirs(maf_dir_out)
    maf_set = set()
    both_set = set()
    maf_exts = []
    for maf_filename in listdir(args.maf_dir):
        maf_parts = maf_filename.split(extsep)
        if "maf" in maf_parts:
            maf_set.add(maf_parts[0])
            maf_exts.append(extsep + extsep.join(maf_parts[1:]))
    for hla_filename in listdir(args.hla_dir):
        hla_basename, hla_ext = path.splitext(hla_filename)
        hla_stripped_basename = hla_basename[:NUM_CHARS_HLA]
        hla_stripped_filename = hla_stripped_basename + hla_ext
        if ".hla" == hla_ext and hla_stripped_basename in maf_set:
            both_set.add(hla_stripped_basename)
            hla_filepath = path.join(args.hla_dir, hla_filename)
            lines = [strip_info(line) for line in open(hla_filepath)]
            hla_filepath_new = path.join(hla_dir_out, 
                    hla_stripped_filename) 
            with open(hla_filepath_new, 'a') as nf:
                output = '\n'.join(lines)
                nf.write(output)
    maf_dirpath = path.abspath(args.maf_dir)
    maf_dirpath_out = path.abspath(maf_dir_out)
    for basename in both_set:
        maf_filepath_original = path.join(maf_dirpath, basename) + maf_exts[0]
        maf_basepath = path.join(maf_dirpath_out, basename)
        maf_filepath = maf_basepath + maf_exts[0]
        maf_filepath_new = maf_basepath + ".maf"
        shutil.copyfile(maf_filepath_original, maf_filepath)
        rename(maf_filepath, maf_filepath_new)

def strip_info(line):
    line = line.strip()
    line = line.split()[0]
    line = line.split('\'')[0]
    return line

if __name__ == "__main__":
    run()

