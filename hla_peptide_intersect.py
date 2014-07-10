#!/usr/bin/env python3
'''
Given a dir_path of allele files (with the first column in each file
being a peptide string) and an output_path for a new CSV file, this
script constructs a DataFrame of the form:

          A0201     A2501     B3501
A0201  1.000000  0.002022  0.014408
A2501  0.002022  1.000000  0.025353
B3501  0.014408  0.025353  1.000000

Each number is the Jaccard set coefficient representing the set 
difference between the peptide strings for the two alleles.

The DataFrame is exported to a CSV file (output_path). Run the
script with --families to group calculations by allele family
(e.g. the set difference between the peptipes for A02* and B35*).

Run the script with --length or --ratio for non-Jaccard calculations.
'''

import pandas as pd
import argparse
from os import path, listdir
import csv

# TODO: Don't double-calculate each comparison (for commutative calculations)
def run(args):
    dir_path = args.dir_path[0]
    calculation = choose_calculation(args)
    files = listdir(dir_path)
    name_to_set = {}
    for file_name in files:
        print('Reading file ' + file_name + '...')
        name = path.splitext(file_name)[0]
        with open(path.join(dir_path, file_name)) as f:
            reader = csv.reader(f, delimiter='\t')
            name_to_set[name] = set([row[0] for row in reader])
    data = {}
    something_to_set = get_something_to_set(name_to_set, 
                                                args.families)
    sorted_items = sorted(something_to_set.items())
    for this_name, this_set in sorted_items:
        print('Doing set comparisons for allele or family ' + this_name)
        data[this_name] = [calculation(this_set, other_set)
                           for other_name, other_set
                           in sorted_items]
    df = pd.DataFrame(data)
    df_columns = df.columns
    df = df.transpose()
    df.columns = df_columns
    df.to_csv(args.output_path[0])
    return df

def choose_calculation(args):
    '''A simple helper method to choose our calculation method.'''
    if args.length:
        return len_intersected
    elif args.ratio:
        return ratio_intersected
    else:
        return jaccard

def get_something_to_set(name_to_set, is_families):
    '''If is_families is true, returns a mapping of allele families
       to sets. Otherwise, returns a mapping of alleles to sets.'''
    if is_families:
        family_to_set = {}
        for allele_name, allele_set in name_to_set.items():
            family_name = allele_name[:3]
            if family_name in family_to_set:
                # TODO: Fix immutable unions
                family_to_set[family_name] = family_to_set[family_name].union(
                        allele_set)
            else:
                family_to_set[family_name] = allele_set
        return family_to_set
    return name_to_set 

def len_intersected(set_a, set_b):
    '''Calculates the length of the set intersection between sets a and b'''
    intersected = set_a.intersection(set_b)
    return len(intersected)

def ratio_intersected(set_a, set_b):
    '''Calculates the length of the set intersection between sets a and b
       over the length of set a'''
    intersected = set_a.intersection(set_b)
    return len(intersected) / len(set_a)

def jaccard(set_a, set_b):
    '''Calculates the Jaccard set difference between sets a and b'''
    intersected = set_a.intersection(set_b)
    unioned = set_a.union(set_b)
    return len(intersected) / len(unioned)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='hla_peptide_intersect.py', 
        usage=__doc__)
    parser.add_argument('dir_path', metavar='dir_path', nargs=1, 
        help='a path to the directory containing allele files')
    parser.add_argument('output_path', metavar='output_path', nargs=1,
        help='the output CSV path')
    parser.add_argument('--families', action='store_true', default=False,
        help='run with --families to group sets by allele families rather')
    parser.add_argument('--length', action='store_true', default=False,
        help='run with --length to count the set intersection length')
    parser.add_argument('--ratio', action='store_true', default=False,
        help='run with --ratio to calculate the set intersection ratio')
    args = parser.parse_args()
    df = run(args)
