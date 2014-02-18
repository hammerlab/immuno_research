#!/usr/bin/env python

"""

Given a directory of refseq protein fasta files and a maf file, generate a pickle file that
maps each maf file to a list of epitopes.

Example:

./maf_to_epitopes.py --refseq-dir data/refseq \ 
					 --maf data/maf/step4_LUSC_Paper_v8.aggregated.tcga.maf2.4.migrated.somatic.maf \
					 epitopes_raw.pkl

 Copyright (c) 2014. Mount Sinai School of Medicine
 
"""


import gzip, logging, argparse, glob, re, pickle

import pandas as pd
import Bio.SeqIO

import common

def refseq_id_to_sequence(refseq_filenames):
	result = {}
	for filename in refseq_filenames:
		logging.info("loading: %s", filename)
		if filename.endswith(".gz"):
			fd = gzip.GzipFile(filename)
		else:
			fd = open(filename)
		gen = Bio.SeqIO.parse(fd, "fasta")
		for record in gen:
			try:
				name = record.id.split("|")[3]
				if name.startswith("NP_"):
					# TODO(odonnt02): handle multiple entries more intelligently than this.
					before_dot = name.split('.')[0]
					result[before_dot] = record
			except IndexError:
				pass
		fd.close()
	logging.info("loaded %d refseq sequences from %d files", len(result), len(refseq_filenames))
	return result

def open_maf(filename):
	logging.info("Opening %s" % filename)
	return pd.read_csv(filename, skiprows=4, sep="\t", low_memory=False)


SINGLE_AMINO_ACID_SUBSTITUTION = re.compile("p.([A-Z])([0-9]+)([A-Z])")
def extract(maf_df, refseq_map, epitope_max_length):
	half_max = int(epitope_max_length / 2)
	filtered = maf_df[["Refseq_prot_Id", "Protein_Change"]].dropna()
	result = []
	for (index, row) in filtered.iterrows():
		ref_seq = refseq_map.get(row.Refseq_prot_Id)
		if ref_seq is None:
			continue
		logging.debug("rfseq match %s", row.Refseq_prot_Id)
		match = SINGLE_AMINO_ACID_SUBSTITUTION.match(row.Protein_Change)
		if match is None:
			continue

		(wild_type, position, mutation) = match.groups()
		position = int(position) - 1

		if wild_type == mutation:
			continue

		if len(ref_seq.seq) <= position:
			logging.warning("Mismatch in protein %s: ref is only %d long, but mutation is at pos %d",
							row.Refseq_prot_Id,
							len(ref_seq.seq),
							position)

			continue
		if ref_seq.seq[position] != wild_type:
			logging.warning("Mismatch in protein %s at pos %d: ref is %s but expected %s",
							row.Refseq_prot_Id,
							position,
							ref_seq.seq[position],
							wild_type)
			continue

		# Make mutation
		full_epitope = ref_seq.seq.tomutable()
		full_epitope[position] = mutation
		trimmed_epitope = full_epitope[max(0, position - half_max) :
									   min(len(full_epitope) - 1, position + half_max)]
		result.append(trimmed_epitope)
	return result

def go():
	parser = argparse.ArgumentParser(usage = __doc__)
	parser.add_argument("--refseq-file", action="append", default=[])
	parser.add_argument('--refseq-dir', action="append", default=[])
	parser.add_argument('--refseq-num', type=int, metavar="X",
		help="Read only first X refseq files (for debugging quickly).")
	parser.add_argument('--maf', action="append", default=[])
	parser.add_argument('--epitope-max-length', type=int, default=30)
	parser.add_argument("out")

	args = parser.parse_args()

	refseq_filenames = list(args.refseq_file)
	for dir in args.refseq_dir:
		refseq_filenames.extend(glob.glob("%s/*" % dir))
	if args.refseq_num:
		refseq_filenames = refseq_filenames[:args.refseq_num]
	refseq_map = refseq_id_to_sequence(refseq_filenames)

	maf_to_epitopes = {}
	for maf_filename in args.maf:
		df = open_maf(maf_filename)
		epitopes = extract(df, refseq_map, args.epitope_max_length)
		logging.info("Matched %d epitopes", len(epitopes))
		maf_to_epitopes[maf_filename] = epitopes

	with open(args.out, 'w') as fd:
		pickle.dump(maf_to_epitopes, fd)
	logging.info("Wrote: %s", args.out)
		




if __name__ == "__main__":
	go()
	





