configfile: 'pangenome_analysis.yaml'

rule all:
    input:
        analysis_dir + "/pangenome/roary/pan_genome_reference.fa"


rule Roary:
    input:
        glob.glob(analysis_dir +"/assemblies/*/annotations/*/*.gff3")
    output:
        outdir = directory(analysis_dir +"/pangenome/roary/"),
        outfile = analysis_dir + "/pangenome/roary/pan_genome_reference.fa"
    threads: 360
    conda: "roary"
    message: "running roary on all of our annotations!"
    log: analysis_dir + "/reports/pangenome/roary.txt"
    shell:
        "roary {input}"

rule FastTree:
    input:
        analysis_dir + "/pangenome/roary/core_gene_alignment.aln"
    output:
        analysis_dir + "/pangenome/roary/tree/core_gene_alignment.newick"
    threads: "360"
    conda: "fasttree"
    message: "creating a newick tree from the output core gene alignment"
    log: analysis_dir + "/reports/pangenome/fasttree.txt"
    shell:
        "FastTree -nt -gtf {input} > {output}"

rule Group2ID:
    input:
    output:
    threads:
    conda: "plasmid_id"
    message:
    log:
    shell:
