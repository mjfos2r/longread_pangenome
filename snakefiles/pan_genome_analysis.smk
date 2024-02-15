configfile: 'configs/pangenome_analysis.yaml'

analysis_dir = config['analysis_dir']

rule all:
    input:
        analysis_dir + "/pangenome/roary/pan_genome_reference.fa",
        analysis_dir + "/pangenome/tree/core_gene_alignment.newick"


rule Roary:
    input:
        glob.glob(analysis_dir +"/assemblies/*/annotations/*/*.gff3")
    output:
        outdir = directory(analysis_dir +"/pangenome/roary/"),
        outfile = analysis_dir + "/pangenome/roary/pan_genome_reference.fa"
    threads: 360
    conda: "roary"
    message: "running roary on all of our annotations!"
    log: analysis_dir + "/logs/pangenome/roary.txt"
    shell:
        "echo 'roary -p 360 -f {output.ourdir} {input} 2>&1 >{log}' "
        "roary -p 360 -f {output.ourdir} {input} 2>&1 >{log}"

rule FastTree:
    input:
        analysis_dir + "/pangenome/roary/core_gene_alignment.aln"
    output:
        analysis_dir + "/pangenome/tree/core_gene_alignment.newick"
    threads: "360"
    conda: "roary"
    message: "creating a newick tree from the output core gene alignment"
    log: analysis_dir + "/logs/pangenome/fasttree.txt"
    shell:
        "FastTree -nt -gtf {input} > {output}"

#rule Group2ID:
#    input:
#    output:
#    threads:
#    conda: "plasmid_id"
#    message:
#    log:
#    shell:
