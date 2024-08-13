# pseudocode/layout for parsing blast hits and determining best hit for each region.


1. parse xml file
2. iterate through records (contigs)
3. iterate through alignments (alignments to reference sequences)
    a. create interval tree for which we will store intervals for each alignment's hsps, then use that tree to determine coverage.
4. iterate through HSPs for each alignment.
    a. add start, end, to our interval tree
    b. determine which homology interval group this hsp is part of
    c. determine percent identity and percent coverage for this hsp
        - query_len = (query_end - query_start+1)
        - query_coverage =  query_len / seq_len * 100
        - query_identity = identities / seq_len * 100

