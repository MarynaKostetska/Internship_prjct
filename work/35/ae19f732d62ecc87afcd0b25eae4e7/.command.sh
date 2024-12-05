#!/bin/bash -ue
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_003663.2&rettype=gb&retmode=text" -O NC_003663.2.gb
