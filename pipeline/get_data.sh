#!/bin/sh

# Download iedb
wget http://www.iedb.org/doc/tcell_compact.zip
unzip tcell_compact.zip

# Download ensembl data

wget ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cds/Homo_sapiens.GRCh37.74.cds.all.fa.gz
gunzip Homo_sapiens.GRCh37.74.cds.all.fa.gz

wget ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/pep/Homo_sapiens.GRCh37.74.pep.all.fa.gz
gunzip Homo_sapiens.GRCh37.74.pep.all.fa.gz

wget ftp://ftp.ensembl.org/pub/release-74/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.74.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh37.74.cdna.all.fa.gz
