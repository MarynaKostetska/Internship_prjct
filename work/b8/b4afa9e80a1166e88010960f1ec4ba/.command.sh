#!/bin/bash -ue
genbank_to -g NC_003663.2.gb -n NC_003663.2.fasta
genbank_to -g NC_003663.2.gb --gff3 NC_003663.2.gff
