import glob

configfile: 'configs/pangenome_analysis.yaml'

analysis_dir = config['analysis_dir']

rule all:
    input:
        analysis_dir + "/pangenome/roary_v2/core_gene_alignment.aln",
        analysis_dir + "/pangenome/tree/v2/fasttree_roary_core_v2.newick",
        analysis_dir + "/pangenome/tree/v2/RAxML_bestTree.roary_v2_gtrgam",
        analysis_dir+"/pangenome/bakta_v2/roary_core_v2.gff3"

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
        "roary -e --mafft -p 360 -f {output.outdir} {input} >{log} 2>&1"

rule FastTree:
    input:
        analysis_dir + "/pangenome/roary_v2/core_gene_alignment.aln"
    output:
        analysis_dir + "/pangenome/tree/v2/fasttree_roary_core_v2.newick"
    threads: 360
    conda: "roary"
    message: "running fasttree on the roary output ==>> {input}"
    log: analysis_dir + "/logs/pangenome/fasttree.txt"
    shell:
        "FastTreeMP -gamma -nt -gtr {input} > {output} 2>{log}"

rule RAxML:
    input:
        analysis_dir + "/pangenome/roary_v2/core_gene_alignment.aln"
    output:
        outdir = directory(analysis_dir + "/pangenome/tree/v2/"),
        outfile = analysis_dir + "/pangenome/tree/v2/RAxML_bestTree.roary_v2_gtrgam"
    threads: 360
    conda: "raxml"
    message: "running RAxML on the roary output ==>> {input}"
    log: analysis_dir + "/logs/pangenome/RAxML_log.txt"
    shell:
        "raxmlHPC -m GTRGAMMA -p 02155 -s {input} -n roary_v2_gtrgam -w {output.outdir} >{log} 2>&1"

#rule Group2ID:
#    input:
#    output:
#    threads:
#    conda: "plasmid_id"
#    message:
#    log:
#    shell:

rule bakta:
    input: analysis_dir + "/pangenome/roary_v2/pan_genome_reference.fa"
    output:
        outdir = directory(analysis_dir+"/pangenome/bakta_v2/"),
        outgff = analysis_dir+"/pangenome/bakta_v2/roary_core_v2.gff3",
    params: bakta_params = config['bakta_params']
    threads: 360
    log: analysis_dir + "/logs/pangenome/bakta_v2.log"
    conda: "bakta"
    message:
        "Running bakta to annotate the pangenome reference at ==>> {input}"
    shell:
        "bakta --db {params.bakta_params[db]} {params.bakta_params[gram]} --threads {threads} {params.bakta_params[opts]} --force --output {output.outdir} --prefix roary_pangene_v2 {input} 2>&1 >{log}"
