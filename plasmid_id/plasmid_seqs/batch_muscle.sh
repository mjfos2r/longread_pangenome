
for file in parA_multi/*.fasta; do
    output_file="alignments/${file%.fasta}.afa"
    muscle -align "$file" -output "$output_file"
done
