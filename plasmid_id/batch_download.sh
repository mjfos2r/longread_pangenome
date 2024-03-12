#!/bin/zsh

#pull in zsh config
source ~/.zshrc
#activate x86 rosetta conda env
condax86
#activate conda env
conda activate entrez

ACCESSION_IDS=$1
OUT=$2

#whereami?
echo "currently in: $PWD"
echo "working with $1 as list of accession ID's"
echo "working with $2 as directory to output fastas into"

echo "lets make sure those output directories exist"
mkdir -p $OUT

echo "okay lets start downloading!"

while read -r line;
do
	echo "downloading file: $line";
	efetch -db nuccore -id $line -format genbank >$OUT/$line.fa;
done < $ACCESSION_IDS

ls -lh $OUT

exit 0
