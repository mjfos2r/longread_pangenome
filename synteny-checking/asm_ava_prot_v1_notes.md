# running big synteny job

1. setup and re-mount results disk.
2. change threads per alignment job from 4 to 2. # and then down to one to maximize the number of parallel jobs :[
3. ran command using new thread specification:
3a. realized that every alignment failed because I didn't have MMSeqs2 installed.
3b. installed MMSeqs2 with AVX2 instructions. __( AND ADDED TO PATH IN .bashrc!!!!! )__
3c. ran the following command:

```{sh}
python3 run_synteny_MP.py --mode ava --input_dir asm_genbanks --output /mnt/disks/results-1/asms_ava_prot_v1 --cores 359 --prog mmseq > asms_ava_prot_v1_run.log
```

4. parsed using command:

```{sh}
python parse_homology_MP_v3.py --dale --version v1 \
--logfile asm_ava_prot \
--annotations_dir asm_genbanks \
--alignments_dir /mnt/disks/results-1/asms_ava_prot_v1/synteny \
--output_dir asms_ava_prot
```

zipping everything for transport to the bucket is as follows:

```{sh}
mkdir -p /mnt/disks/results-1/asms_ava_prot_v1/tarballs && \
find /mnt/disks/results-1/asms_ava_prot_v1/synteny/ -mindepth 1 -maxdepth 1 -type d -print0 | \
parallel -0 -j 4 'tar -czf "/mnt/disks/results-1/asms_ava_prot_v1/tarballs/$(basename {}).tar.gz" -C "$(dirname {})" "$(basename {})"'
```

and then parallel zipping all the tarballs with pigz is as follows:

```{sh}
find tarballs -type f -print0 | parallel -0 -N1000 -P$(nproc) tar cf - --no-recursion {} | pigz > asms_ava_gene_homology_v1.tar.gz
```

then move that tarball to the bucket! it is stored as a "composite object" which means that we'll cross that bridge when we come to it!