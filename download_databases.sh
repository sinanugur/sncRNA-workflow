#!/bin/bash




if [ ! -f databases/tRNA.gff3 ]; then

latest=$(curl -s4 anonymous@ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/  | grep -Eo 'release_[0-9]+' | sort -nrt _ -k2,2 | head -n 1)
release=$(echo $latest | cut -d_ -f2)


rm databases/tRNA.gff3
rm databases/gencode.gff3.gz


wget "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/$latest/gencode.v$release.tRNAs.gff3.gz" -O databases/gencode.v$release.tRNAs.gff3.gz

gzcat databases/gencode.v$release.tRNAs.gff3.gz > databases/tRNA.gff3


wget  "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/$latest/gencode.v$release.annotation.gff3.gz" -O databases/gencode.v$release.annotation.gff3.gz

mv databases/gencode.v$release.annotation.gff3.gz databases/gencode.gff3.gz

echo $release > databases/dbversion.txt

fi


if [ ! -f databases/hg38.fa.gz ]; then

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
mv hg38.fa.gz databases/
fi
