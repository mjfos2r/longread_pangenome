import glob

configfile: 'configs/pangenome_analysis.yaml'

analysis_dir = config['analysis_dir']

rule all:
    input:
        analysis_dir + "/pangenome/roary_v2/core_gene_alignment.aln",
        analysis_dir + "/pangenome/tree/fasttree_roary_core_v2.newick",
        analysis_dir + "/pangenome/tree/RAxML_bestTree_roary_v2_gtrgam"

rule Roary:
    input:
        glob.glob(analysis_dir +"/assemblies/*/annotation/*/*.gff3")
    output:
        outdir = directory(analysis_dir +"/pangenome/roary_v2/"),
        outfile = analysis_dir + "/pangenome/roary_v2/core_gene_alignment.aln"
    threads: 360
    message: "running roary on all of our annotations! : \n roary -p 360 -f {output.outdir} {input} 2>&1 >{log}"
    log: analysis_dir + "/logs/pangenome/roary_log.txt"
    conda: "roary"
    shell:
        "roary -e -p 360 -f {output.outdir} {input} >{log} 2>&1"

rule FastTree:
    input:
        analysis_dir + "/pangenome/roary_v2/core_gene_alignment.aln"
    output:
        analysis_dir + "/pangenome/tree/fasttree_roary_core_v2.newick"
    threads: 360
    conda: "roary"
    message: "running fasttree on the roary output ==>> {input}"
    log: analysis_dir + "/logs/pangenome/fasttree.txt"
    shell:
        "FastTree -gamma -nt -gtr {input} > {output} 2>{log}"

rule RAxML:
    input:
        analysis_dir + "/pangenome/roary_v2/core_gene_alignment.aln"
    output:
        analysis_dir + "/pangenome/tree/RAxML_bestTree_roary_v2_gtrgam"
    threads: 360
    conda: "raxml"
    message: "running RAxML on the roary output ==>> {input}"
    log: analysis_dir + "/logs/pangenome/RAxML_log.txt"
    shell:
        "raxmlHPC -m GTRGAMMA -p 02155 -s {input} -n roary_v2_gtrgam >{log} 2>&1"

#rule Group2ID:
#    input:
#    output:
#    threads:
#    conda: "plasmid_id"
#    message:
#    log:
#    shell:
