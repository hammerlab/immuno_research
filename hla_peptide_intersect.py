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

Run the script with --myfamily to compare each allele with its own
allele family (minus itself).

Run the script with --length or --ratio for non-Jaccard calculations.
'''

import pandas as pd
import argparse
from os import path, listdir
import csv
from collections import defaultdict

def run(args):
    dir_path = args.dir_path[0]
    calculation = choose_calculation(args)
    files = listdir(dir_path)
    allele_to_set = {}
    for file_name in files:
        print('Reading file ' + file_name + '...')
        name = path.splitext(file_name)[0]
        with open(path.join(dir_path, file_name)) as f:
            reader = csv.reader(f, delimiter='\t')
            allele_to_set[name] = set([row[0] for row in reader])
    data = {}
    if args.myfamily:
        allele_to_family_set = get_allele_to_family_set(allele_to_set)
        for allele_name, family_set in \
                sorted(allele_to_family_set.items()):
            print('Doing set comparisons for allele with family ' + 
                allele_name)
            data[allele_name] = [calculation(allele_to_set[allele_name], 
                                             family_set)]
        df = pd.DataFrame(data)
        df = df.transpose()
        df.to_csv(args.output_path[0])
        return df
    something_to_set = get_family_to_set(allele_to_set) if \
            args.families else allele_to_set
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

def get_family_to_set(allele_to_set):
    '''Returns a mapping of allele families to sets.'''
    family_to_set = defaultdict(set)
    for allele_name, allele_set in allele_to_set.items():
        family_name = allele_name[:3]
        if family_name in family_to_set:
            family_to_set[family_name].update(allele_set)
    return family_to_set

def get_family_to_names(allele_to_set):
    '''Return a mapping of allele familes to allele names in that 
       family.'''
    family_to_names = defaultdict(set)
    for allele_name in allele_to_set.keys():
        family_name = allele_name[:3]
        family_to_names[family_name].add(allele_name)
    return family_to_names

def get_allele_to_family_set(allele_to_set):
    '''Returns a mapping of allele to its family's set (without it).'''
    allele_to_family_set = defaultdict(set)
    family_to_names = get_family_to_names(allele_to_set)
    for allele_name, allele_set in allele_to_set.items():
        family_name = allele_name[:3]
        sibling_names = family_to_names[family_name].copy()
        sibling_names.remove(allele_name)
        for sibling_name in sibling_names:
            allele_to_family_set[allele_name].update(
                    allele_to_set[sibling_name])
    return allele_to_family_set

def len_intersected(set_a, set_b):
    '''Calculates the length of the set intersection between sets a
       and b'''
    intersected = set_a.intersection(set_b)
    return len(intersected)

def ratio_intersected(set_a, set_b):
    '''Calculates the length of the set intersection between sets a 
       and b over the length of set a'''
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
        help='run with --families to group sets by allele families')
    parser.add_argument('--length', action='store_true', default=False,
        help='run with --length to count the set intersection length')
    parser.add_argument('--ratio', action='store_true', default=False,
        help='run with --ratio to calculate the set intersection ratio')
    parser.add_argument('--myfamily', action='store_true', default=False,
        help='run with --myfamily to compare alleles with their families')
    args = parser.parse_args()
    df = run(args)
