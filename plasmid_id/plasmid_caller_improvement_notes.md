# Finally Fix Plasmid Caller

1.  Get assembly.
2.  Get contigs from assembly. (parse gbff)
4.  Get genes on contig. BioPython
5.  BLAST vs PF32_db
6.  BLAST vs WP_db
7.  Parse BLAST hits
8.  Determine best call
9.  Take best-call, (use WP best ident?)
10. Get genes from best call
11. Compare gene presence/order for contig against best-call
12. Determine presence/absence of genes on contig
13. Determine order of contig vs best-call
14. Compile pf32 best-call, wp best-call, gene order best-call.
15. Final call determined from gene order. 


# Sanity Check of Caller

>> this is basically just 10 - 13 from above.
>> we can use Rachel's tool to generate heatmaps for gene presence across isolates.
>> Take ref plasmid, get genes, put in order.
>> Take contig, get genes, compare to ref.
>> Plot, determine how scuffed it is.