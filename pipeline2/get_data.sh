#!/bin/sh

mkdir -p data/refseq
mkdir -p data/maf

# MAF files of interest
wget --directory-prefix=data/maf -nc \
	https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/lusc/gsc/broad.mit.edu/illuminaga_dnaseq/mutations/broad.mit.edu_LUSC.IlluminaGA_DNASeq.Level_2.100.1.0/step4_LUSC_Paper_v8.aggregated.tcga.maf2.4.migrated.somatic.maf 

# Full proteome
wget --directory-prefix=data/refseq -nc \
	ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/*.faa.gz
