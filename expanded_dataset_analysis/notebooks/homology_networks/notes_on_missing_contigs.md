# Reviewing the AVA alignment that I generated using the below mummer command:

Where `$1` is the multifasta, `$2` is the output prefix, and `$3` is the minimum length.
```bash
docker run -v $(pwd):/data -it mjfos2r/mummer4 -run /data/assemblies/all_contigs_v5.fna /data/output/alignments/ava/new_ava/nucl_v5 250
```

## Okay so let's look at our input 
and see how many contigs we're running an alignment of.
```bash
cat ../../../../assemblies/all_contigs_v5.fna | grep ">" | wc -l
4080
```
## okay so we have 4080 contigs. 
but when we look at unique query ids in our mummer output, we only get ~2698 nodes. 
```bash
parse_mummer_output.py output/alignments/ava/new_ava/nucl_v5.coords.tab

Network Statistics:
Total edges: 258997
Unique nodes: 2698
Forward alignments: 0
```
## But if we review the details of our contigs, 
running the following to filter contigs by a 250bp cutoff
```python
from Bio import SeqIO
# lists for our short and long contigs
under_250 = []
over_250 = []
# Read and filter
for record in SeqIO.parse("../../assemblies/all_contigs_v5.fna", "fasta"):
    if len(record.seq) < 250:
        under_250.append(record.id)
    else:
        over_250.append(record.id)
# summary
print(f"Contigs under 250bp: {len(under_250)}")
print(f"Contigs over 250bp: {len(over_250)}")
```
which returns:
```bash
Contigs under 250bp: 1221
Contigs over 250bp: 2859
```
## Okey dokey so far this tracks.
4080 contigs total
2859 contigs >= 250bp
2698 contigs in homology network.
## Let's see if we can bump that up by dropping the minimum alignment length down to 200 bp and then re-running the above.
***
## ATTENTION: nucl_v5 HAS BEEN RENAMED TO nucl_v5_250 to keep things groovy.
- Also, reviewing the documentation for nucmer, setting the minimum length is unnecessary as it is default to 20bp. We're just gonna strip `$3` and the `-l` flag from the script and let'er'rip.
***
## RE-RUN
ok let's try again
### running with min_alignment_length set to 200.
```bash
docker run -v $(pwd):/data -it mjfos2r/mummer4 -run /data/assemblies/all_contigs_v5.fna /data/output/alignments/ava/new_ava/nucl_v5_20
```
### network re-gen
```bash
parse_mummer_output.py output/alignments/ava/new_ava/nucl_v5_20.coords.tab

```
which gives:
```bash
Network Statistics:
Total edges: 866725
Unique nodes: 3960
Forward alignments: 0
Reverse alignments: 866725
Average alignments per edge: 8.34
```
### Ok groovy. Let's check the actual output network file and make sure for ourselves.
```bash
awk '{print $1}' nucl_v5_20_homology_network.tsv | uniq | wc -l
3921
```
HM SOMETHING DOESNT ADD UP HERE.

It appears that .ipynb_checkpoints polluted my multifasta by including a checkpoint assembly:
`ESI26H-checkpoint__contig##`
unacceptable.
I will strip it from the multifasta:
```bash
awk '/^>/ {p=!/-checkpoint/} p' all_contigs_v5.fna > all_contigs_v5_fixed.fna
```
ATTN: THIS PROBABLY NEEDS TO JUST BE GENERATED FROM THE BEGINNING USING THE FIXED FILE!!
I will strip it from the coords file at once and regenerate the network.
```bash
awk '!/-checkpoint/' OFS='\t' nucl_v5_20.coords.tab > fixed_nucl_v5_20.coords.tab
awk '!/-checkpoint/' OFS='\t' nucl_v5_20.coords > fixed_nucl_v5_20.coords
```
re-gen gives:
```bash
Network Statistics:
Total edges: 842549
Unique nodes: 3894
Forward alignments: 0
Reverse alignments: 842549
Average alignments per edge: 8.45
```
I'm so damn confused.
```bash
awk -F'\t' '{print $1}' homology_network.tsv | sort | uniq | wc -l
3809
```
What 85 nodes is this script hallucinating?

*** 
Ok let's rewrite the parser using this handy info from the mummer docs. 
```bash
"When run with the -B option, output format will consist of 21 tab-delimited columns. These are as follows: 
[1] query sequence ID 
[2] date of alignment 
[3] length of query sequence 
[4] alignment type 
[5] reference file 
[6] reference sequence ID 
[7] start of alignment in the query 
[8] end of alignment in the query 
[9] start of alignment in the reference 
[10] end of alignment in the reference 
[11] percent identity 
[12] percent similarity 
[13] length of alignment in the query 
[14] 0 for compatibility 
[15] 0 for compatibility 
[16] NULL for compatibility 
[17] 0 for compatibility 
[18] strand of the query 
[19] length of the reference sequence 
[20] 0 for compatibility 
[21] and 0 for compatibility."
```

after substantially rewriting the parsing. We get the following: 
```bash
python ~/longread_pangenome/expanded_dataset_analysis/scripts/parse_mummer_output.py fixed_nucl_v5_20.coords.tab

Network Statistics:
Total edges: 466494
Unique nodes: 3908
Average alignments per edge: 8.04

Alignment Length Statistics:
Total alignment length: 5,073,393,965
Mean alignment length: 10875.58
Median alignment length: 769.00
```

which like, ok. internal cutoffs are for 250bp and 90% ident. Edges don't care about fwd or rev alignments.
everything's piled up and squished into a single edge to give us the aln len.
Lets throw this into a matrix and then plot our graph now.

## Matrix Maker
ok let's make this thing.
```bash
python ~/longread_pangenome/expanded_dataset_analysis/scripts/make_homology_matrix.py fixed_nucl_v5_20_homology_network.tsv ~/longread_pangenome/expanded_dataset_analysis/assemblies/all_contigs_v5_fixed.fna

Matrix Statistics:
Matrix dimensions: 3908 x 3908
Sum of diagonal elements (total contig lengths): 115,614,563
Number of non-zero off-diagonal elements: 932988
Average non-zero alignment length: 10875.58
```
*_groovy_*

okey now let's feed this into the graph plotter notebook!!!

## GRAPH PLOTTIN.

***
## OKAY SO I NEED TO GO BACK AND FILTER AT 1000BP SINCE THAT'S THE MIN CUTOFF FOR OUR CLASSIFICATION METHODS.
### OR DO I ??
***
Nah we're gonna let it rip at 250bp. It looks *very* weird at 1000bp cutoff. I'll need to review the old code.

***

**ATTENTION: new_ava has been renamed to new_ava_v1!**

***
# Okay round 3 let's get this functional.
Actually generating the final plot.
## MUMmer AVA run
Starting all the way from the multiple alignment mummer run.

Command:
```bash
docker run -v $(pwd):/data -it mjfos2r/mummer4 \
-run \
/data/assemblies/all_contigs_v5_fixed.fna \
/data/output/alignments/ava/new_ava_v3/fixed_nucl_v5
```
Output:
```bash
real	21m56.244s
user	607m7.419s
sys	0m11.217s
WARNING: the options -bdHcloTw are all ignored for btab
```
***
## Parse alignments into network
ok let's now parse this thing and generate our network table. 
Length and identity cutoffs changed in script to: 750bp, 85%

Command:
```bash
python scripts/parse_mummer_output.py output/alignments/ava/new_ava_v3/fixed_nucl_v5.coords.tab output/homology_networks/nucl_v5_v3/
```
Output:
```bash
Network Statistics:
Total edges: 270983
Unique nodes: 2829
Average alignments per edge: 12.25

Alignment Length Statistics:
Total alignment length: 5,113,864,881
Mean alignment length: 18871.53
Median alignment length: 2770.00
```

## Parse network table into matrix
ok let's parse that network into a matrix.

Command:
```bash
python scripts/make_homology_matrix.py output/homology_networks/nucl_v5_v3/fixed_nucl_v5_homology_network.tsv assemblies/all_contigs_v5_fixed.fna output/homology_networks/nucl_v5_v3/
```
Output:
```bash
Matrix Statistics:
Matrix dimensions: 2829 x 2829
Sum of diagonal elements (total contig lengths): 115,433,382
Number of non-zero off-diagonal elements: 541966
Average non-zero alignment length: 18871.53
```

## Ok now let's run our generation script. {{TODO: convert notebook for graph generation into script!}}
For context, the previous dataset was filtered at 1000bp and had the following number of nodes:
```bash
Total vertices: 1100
Aligned labels: 1100
Aligned groups: 1100
Unique groups: {
    'lp28-7',  'cp32-6',     'lp36',     'lp28-2', 
    'lp17',    'cp9',        'lp28-9',   'lp54', 
    'cp32-13', 'cp32-1',     'cp32-12',  'lp21', 
    'cp32-3',  'chromosome', 'lp56',     'cp9-3', 
    'cp32-7',  'lp28-6',     'lp28-1',   'lp28-4', 
    'lp38',    'cp32-9',     'lp28-11',  'cp32-10', 
    'cp32-4',  'cp32-3+10',  'cp32-11',  'cp32-5', 
    'lp5',     'lp28-5',     'lp21-cp9', 'lp28-3', 
    'cp32-8',  'lp25',       'lp28-8',   'cp26',
}
```
Ok still running into issues with the weights being >1.0, this is because the alignment lengths are all wonky because I'm using `--nosimplify`.

I can apparently omit that flag and it will keep alignments for multiple queries but simplify alignments for the same query. 

Ok let's go, from the top!
-- i am going to make a script to auto run everything in line ;) that's done, going to run this manually and then test the script. --

## RERUN MANUALLY BEFORE SSCRIPT TESTING.
```bash
docker run -v $(pwd):/data -it mjfos2r/mummer4 -run /data/assemblies/all_contigs_v5_fixed.fna /data/output/alignments/ava/new_ava_v3/fixed_nucl_v5 #let's try this again ;)

real	18m6.910s
user	515m51.076s
sys	0m18.910s
```

```bash
python scripts/parse_mummer_output.py output/alignments/ava/new_ava_v3/fixed_nucl_v5.coords.tab output/homology_networks/nucl_v5_v3/

Network Statistics:
Total edges: 271210
Unique nodes: 2829
Average alignments per edge: 11.03

Alignment Length Statistics:
Total alignment length: 5,116,242,672
Mean alignment length: 18864.51
Median alignment length: 2774.00
```

```bash
python scripts/make_homology_matrix.py output/homology_networks/nucl_v5_v3/fixed_nucl_v5_homology_network.tsv assemblies/all_contigs_v5_fixed.fna output/homology_networks/nucl_v5_v3/

Matrix Statistics:
Matrix dimensions: 2829 x 2829
Sum of diagonal elements (total contig lengths): 115,433,382
Number of non-zero off-diagonal elements: 542420
Average non-zero alignment length: 18864.51
```
Still getting that weird weight issue.

Gotta be something funky with some alignment somewhere.
Probably not taking Q vs R and R vs Q into account properly.

## Testing pipeline script! 
Ok let's go.
Running with the following: 
```bash
scripts/homology_network_graph_pipeline.sh -i assemblies/all_contigs_v5_fixed.fna -o output -p nucl_v5 -v v4 -m output/genotyping/replicons/calls_v10/dataset_v5_best_replicon_hits.csv
Starting pipeline with:
  Input file: assemblies/all_contigs_v5_fixed.fna
  Output directory: output
  Prefix: nucl_v5
  Version: v4
  ID mapping file: output/genotyping/replicons/calls_v10/dataset_v5_best_replicon_hits.csv
Running pipeline steps...
Parsing MUMmer output...

Network Statistics:
Total edges: 271210
Unique nodes: 2829
Average alignments per edge: 11.03

Alignment Length Statistics:
Total alignment length: 5,116,242,672
Mean alignment length: 18864.51
Median alignment length: 2774.00
Generating homology matrix...

Matrix Statistics:
Matrix dimensions: 2829 x 2829
Sum of diagonal elements (total contig lengths): 115,433,382
Number of non-zero off-diagonal elements: 542420
Average non-zero alignment length: 18864.51
Generating 3D visualizations...
Processing matrix rows: 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 2829/2829 [01:06<00:00, 42.58it/s]

Edge Statistics:
Total edges: 207691
Unique source nodes: 1845
Unique target nodes: 1847
Weight range: 0.002 - 1.617
Saved edges to output/homology_networks/v4/nucl_v5_matrix_edges.csv
Sample of idmap:
  assembly_id   contig_id  contig_len                    plasmid_id plasmid_name  ... covered_intervals     query_intervals  subject_hit_coords  query_coverage_percent  call_method
0       B331P   contig_1       903654        B500_chromosome_ParA_2   chromosome  ...        [(1, 364)]  [(368951, 370042)]          [(1, 364)]                0.120732         pf32
1       B331P  contig_10         8714                gb|CP017210.1|          cp9  ...       [(1, 8714)]         [(1, 8714)]         [(1, 8714)]               99.988524           wp
2       B331P  contig_11        15594  RS00875_MM1_plsm_lp17_ParA_X         lp17  ...        [(1, 242)]    [(12036, 12761)]          [(1, 242)]                4.649224         pf32
3       B331P  contig_12        24722   RS00040_ZS7_ZS7_lp25_ParA_X         lp25  ...        [(1, 252)]    [(11900, 12655)]          [(1, 252)]                3.053960         pf32
4       B331P  contig_13        27641         H28_B31_lp28-3_ParA_X       lp28-3  ...        [(1, 251)]    [(16701, 17453)]          [(1, 251)]                2.720596         pf32

[5 rows x 16 columns]

Sample of vertex_order:
['URI101H__contig_19', 'GCF_040790765.1_ASM4079076v1_genomic__NZ_CP161108.1', 'URI118H__contig_1', 'URI42H__contig_27', 'URI112H__contig_18', 'UCT96H__contig_17', 'UWI263P__contig_14', 'UWI247P__contig_14', 'UNY203P__contig_14', 'URI120H__contig_13']
Fragment found: UCT35H__contig_17 (assigned as unknown)

Total vertices: 1861
Aligned labels: 1861
Aligned groups: 1861
Unique groups: {'cp32-12', 'lp38', 'cp32-6', 'lp36', 'cp32-7', 'unknown', 'cp32-11', 'lp28-6', 'lp17', 'cp26', 'chromosome', 'cp32-1', 'cp32-9', 'cp32-4', 'lp25', 'lp28-2', 'lp32-3', 'lp21-cp9', 'cp32-5', 'lp28-8', 'lp28-4', 'cp9', 'lp28-5', 'lp5', 'lp21', 'lp54', 'cp32-8', 'lp28-11', 'lp28-3', 'lp56', 'cp32-10', 'lp28-7', 'lp28-1', 'cp32-13', 'cp9-3', 'lp28-9', 'cp32-3'}
PNG written to igraph_asm_ava_homology_nucl_kk3d.png
HTML written to output/homology_networks/v4/plots/igraph_asm_ava_homology_nucl_kk3d.html
PNG written to igraph_asm_ava_homology_nucl_kk3d_no_edges.png
Network data saved to output/homology_networks/v4/plots/network.json
Pipeline completed successfully!
Results can be found in: output
```

***
## RERUNNING SCRIPT YET AGAIN BUT THE LAST TIME THIS TIME I PROMISE
Okay I fixed the interval merging step in parse_MUMmer. I
I also added the correct color mapping from Rachel.
I also added the replicon names to the edges table for use in gephi/cytoscape.

Reran the script using the alignments generated for 4_1. 

Output is as follows:
```bash
scripts/homology_network_graph_pipeline.sh -i assemblies/all_contigs_v5_fixed.fna -o output -p nucl_v5 -v v4_2 -m output/genotyping/replicons/calls_v10/dataset_v5_best_replicon_hits.csv
Starting pipeline with:
  Input file: assemblies/all_contigs_v5_fixed.fna
  Output directory: output
  Prefix: nucl_v5
  Version: v4_2
  ID mapping file: output/genotyping/replicons/calls_v10/dataset_v5_best_replicon_hits.csv
  Minimum Length: 1000
  Minimum Weight: 0.5
Running pipeline steps...
Parsing MUMmer output...

Network Statistics:
Total edges: 1014156
Unique nodes: 3911
Average alignments per edge: 6.80

Alignment Length Statistics:
Total alignment length: 8,949,770,814
Mean alignment length: 8824.85
Median alignment length: 871.00
Generating homology matrix...

Matrix Statistics:
Matrix dimensions: 3911 x 3911
Sum of diagonal elements (total contig lengths): 115,624,152
Number of non-zero off-diagonal elements: 2028312
Average non-zero alignment length: 8824.85
Generating 3D visualizations...
/home/mf019/longread_pangenome/expanded_dataset_analysis/scripts/draw_3D_graph.py:24: SettingWithCopyWarning:
A value is trying to be set on a copy of a slice from a DataFrame.
Try using .loc[row_indexer,col_indexer] = value instead

See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
  fixed_map['name'] = fixed_map.apply(lambda x: x['assembly_id'] + '__' + x['contig_id'], axis=1)
Processing matrix rows: 100%|███████████████████████████████████████████████████████| 3911/3911 [02:08<00:00, 30.42it/s]

Edge Statistics:
Total edges: 170689
Unique source nodes: 1621
Unique target nodes: 1618
Weight range: 0.500 - 1.673
Saved edges to output/homology_networks/v4_2/nucl_v5_matrix_edges.csv
Sample of idmap:
  assembly_id  contig_id  contig_len  ... subject_hit_coords query_coverage_percent call_method
0       B331P   contig_1      903654  ...         [(1, 364)]               0.120732        pf32
1       B331P  contig_10        8714  ...        [(1, 8714)]              99.988524          wp
2       B331P  contig_11       15594  ...         [(1, 242)]               4.649224        pf32
3       B331P  contig_12       24722  ...         [(1, 252)]               3.053960        pf32
4       B331P  contig_13       27641  ...         [(1, 251)]               2.720596        pf32

[5 rows x 16 columns]

Sample of vertex_order:
['URI46H__contig_15', 'URI118H__contig_2', 'GCF_040790735.1_ASM4079073v1_genomic__NZ_CP161056.1', 'GCF_040790775.1_ASM4079077v1_genomic__NZ_CP161189.1', 'URI48H__contig_3', 'UNY1038P__contig_3', 'UWI263P__contig_10', 'GCF_002151505.1_ASM215150v1_genomic__NZ_CP019763.1', 'UNY1090P__contig_5', 'URI40H__contig_4']

Total vertices: 1673
Aligned labels: 1673
Aligned groups: 1673
Unique groups: {'lp28-11', 'lp28-1', 'cp32-8', 'cp32-12', 'cp9-3', 'chromosome', 'cp32-13', 'lp5', 'cp9', 'lp25', 'lp36', 'lp56', 'cp32-6', 'lp28-3', 'cp32-7', 'lp28-5', 'lp21-cp9', 'cp32-9', 'cp26', 'cp32-5', 'lp17', 'cp32-3', 'lp28-4', 'cp32-1', 'lp38', 'lp54', 'lp21', 'lp28-8', 'lp28-7', 'lp28-2', 'cp32-4', 'lp28-6', 'cp32-10', 'lp28-9', 'cp32-11'}
PNG written to igraph_asm_ava_homology_nucl_kk3d.png
HTML written to output/homology_networks/v4_2/plots/igraph_asm_ava_homology_nucl_kk3d.html
PNG written to igraph_asm_ava_homology_nucl_kk3d_no_edges.png
Network data saved to output/homology_networks/v4_2/plots/network.json
Pipeline completed successfully!
Results can be found in: output
```

Okay so I just realized that the duplicate clique thing is strictly because of duplicated nodes. This should be fairly straightforward to repair..

Also, my contig id mapping is incorrect, I'm using ben's generated table from now on since his is accurate AND I don't care about contigs < 1000bp.

## STILL BROKEN. RUNNING PGGB NOW. IF ALL ELSE FAILS I'LL USE THE V4_2 NETWORK I'VE ALREADY GENERATED...

```bash
pggb -x auto \
-in /data/assemblies/all_contigs_v5_fixed.fna \
-o /data/output/homology_network/pggb_v1 \
-t 30 \
-p 90 \
-s 1000 \
-V 'ref:1000'
```
*** 
## That did not generate anything usable, I am redoing the network from scratch (parsing MUMmer output that is)

okay so I wrote a solution using intervaltree.

this is running and might take a while so once it's finished, we'll generate our nodes and edges table,
then re-generate the network!