"""

 Copyright (c) 2014. Mount Sinai School of Medicine
 
"""
from pipeline import ImmunoPipeline
import argparse
from immunogenicity import ImmunogenicityRFModel
from binding import IEDBMHCBinding
from cleavage import ProteasomalCleavage
import reduced_alphabet
from Bio import SeqIO
import pandas as pd

def get_epitopes_from_fasta(fasta_file):
  epitope_data = SeqIO.parse(fasta_file, 'fasta')
  return pd.Series(list(set([e.seq for e in epitope_data])))


def create_pipeline():
  pipeline = ImmunoPipeline()
  return pipeline

def add_scoring(pipeline, alleles):
  pipeline.add_scorer(IEDBMHCBinding(name='mhc', alleles=alleles))
  pipeline.add_scorer(ProteasomalCleavage(name='protocleave'))
  #pipeline.add_scorer(SelfScorer(name='selfcheck'))
  pipeline.add_scorer(ImmunogenicityRFModel(name='default RF'))
  pipeline.add_scorer(ImmunogenicityRFModel(name = 'murphy10 RF', reduced_alphabet = reduced_alphabet.murphy10))

  return pipeline

DEFAULT_ALLELE = 'HLA-A*24:92'

if __name__ == '__main__':

  parser = argparse.ArgumentParser()
  parser = argparse.ArgumentParser()
  parser.add_argument("--input", help="input file to process")
  parser.add_argument("--allele_file", help="list of alleles")
  parser.add_argument("--output", help="output file for dataframes", required=True)


  args = parser.parse_args()
  if args.input.endswith(".vcf"):
    epitope_data = pipeline.add(Variant2Epitope)
  elif args.input.endswith(".maf"):
    epitope_data = pipeline.add(MafEpitope)
  elif args.input.endswith(".fasta") or args.input.endswith(".fa"):
    epitope_data =  get_epitopes_from_fasta(args.input)

  if args.allele_file:
    alleles = [l.strip() for l in open(allele_file)]
  else:
    alleles = [DEFAULT_ALLELE]
  pipeline = ImmunoPipeline()
  add_scoring(pipeline, alleles)
  data = pipeline.score(epitope_data)

  data.to_csv(args.output, index=False)
