#!/usr/bin/env python

"""

Given a fasta file of a proteome, write out a pickle file with all the epitopes of the given length.

./generate_self_epitopes.py --size 9 --size 10 --size 11 --size 12 \
	--fasta data/HUMAN.fasta.gz \
  self_epitopes.csv

With --size 10 on the human human proteome, this takes ~2 minutes and writes out a 144 mb file with 
~11.5 million self epitopes.

 Copyright (c) 2014. Mount Sinai School of Medicine
 
"""

import gzip, logging, argparse, pandas, collections
import Bio.SeqIO
import common

def fasta_to_epitopes(filenames, sizes):
  epitopes = collections.Counter()
  for filename in filenames:
	logging.info("loading: %s", filename)
	if filename.endswith(".gz"):
		fd = gzip.GzipFile(filename)
	else:
		fd = open(filename)
	gen = Bio.SeqIO.parse(fd, "fasta")
	for record in gen:
		seq = str(record.seq)
		for size in sizes:
			for i in range(len(record.seq) - size + 1):
				epitopes[seq[i:i+size]] += 1
	fd.close()
	return pandas.DataFrame(epitopes.most_common(), columns=["Epitope", "Num Occurrences"])

def go():
  parser = argparse.ArgumentParser(usage = __doc__)
  parser.add_argument("--size", required=True, type=int, action="append", default=[])
  parser.add_argument('--fasta', required=True, action="append", default=[])
  parser.add_argument("out")

  args = parser.parse_args()

  epitopes_df = fasta_to_epitopes(args.fasta, args.size)
  epitopes_df.to_csv(args.out, index = False)
  logging.info("Wrote %d unmutated epitopes to: %s", len(epitopes_df), args.out)
    
if __name__ == "__main__":
  go()
  
