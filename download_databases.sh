#!/bin/bash


latest=$(curl -s4 anonymous@ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/  | grep -Eo 'release_[0-9]+' | sort -nrt _ -k2,2 | head -n 1)
release=$(echo $latest | cut -d_ -f2)





wget "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/$latest/gencode.v$release.tRNAs.gff3.gz" -O databases/gencode.v$release.tRNAs.gff3.gz

gzcat databases/gencode.v$release.tRNAs.gff3.gz > databases/tRNA.gff3


wget  "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/$latest/gencode.v$release.annotation.gff3.gz" -O databases/gencode.v$release.annotation.gff3.gz


#gzcat databases/gencode.v$release.annotation.gff3.gz > databases/gencode.gff3

mv databases/gencode.v$release.annotation.gff3.gz databases/gencode.gff3.gz


