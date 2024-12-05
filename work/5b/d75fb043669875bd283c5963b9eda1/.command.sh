#!/bin/bash -ue
awk '!/^>/ {print}' NC_003663.2.fasta > "NC_003663.2_sequence.fasta"
