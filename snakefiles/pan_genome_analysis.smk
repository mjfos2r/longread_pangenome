import glob

configfile: 'configs/pangenome_analysis.yaml'

analysis_dir = config['analysis_dir']

rule all:
    input:
        analysis_dir + "/pangenome/roary_v3/core_gene_alignment.aln",
        analysis_dir + "/pangenome/tree/v3/fasttree_roary_core_v3.newick",
        analysis_dir + "/pangenome/tree/v3/RAxML_bestTree.roary_v3_gtrgam",
        analysis_dir+"/pangenome/bakta_v3/roary_core_v3.gff3"

rule Roary:
    input:
        glob.glob(analysis_dir +"/assemblies/*/annotation/*/*.gff3")
    output:
        outdir = directory(analysis_dir +"/pangenome/roary_v3/"),
        outfile = analysis_dir + "/pangenome/roary_v3/core_gene_alignment.aln"
    threads: 15
    message: "running roary on all of our annotations! : \n roary -p 360 -f {output.outdir} {input} 2>&1 >{log}"
    log: analysis_dir + "/logs/pangenome/roary_log.txt"
    singularity: 'docker://sangerpathogens/roary'
    shell:
        "roary -e --mafft -p 360 -f {output.outdir} {input} >{log} 2>&1"

rule FastTree:
    input:
        analysis_dir + "/pangenome/roary_v3/core_gene_alignment.aln"
    output:
        analysis_dir + "/pangenome/tree/v3/fasttree_roary_core_v3.newick"
    threads: 15
    conda: "fasttree"
    message: "running fasttree on the roary output ==>> {input}"
    log: analysis_dir + "/logs/pangenome/fasttree.txt"
    shell:
        "FastTreeMP -gamma -nt -gtr {input} > {output} 2>{log}"

rule RAxML:
    input:
        analysis_dir + "/pangenome/roary_v3/core_gene_alignment.aln"
    output:
        outdir = directory(analysis_dir + "/pangenome/tree/v3/"),
        outfile = analysis_dir + "/pangenome/tree/v3/RAxML_bestTree.roary_v3_gtrgam"
    threads: 15
    conda: "raxml"
    message: "running RAxML on the roary output ==>> {input}"
    log: analysis_dir + "/logs/pangenome/RAxML_log.txt"
    shell:
        "raxmlHPC -m GTRGAMMA -p 02155 -s {input} -n roary_v3_gtrgam -w {output.outdir} >{log} 2>&1"

#rule Group2ID:
#    input:
#    output:
#    threads:
#    conda: "plasmid_id"
#    message:
#    log:
#    shell:

rule bakta:
    input: analysis_dir + "/pangenome/roary_v3/pan_genome_reference.fa"
    output:
        outdir = directory(analysis_dir+"/pangenome/bakta_v3/"),
        outgff = analysis_dir+"/pangenome/bakta_v3/roary_core_v3.gff3",
    params: bakta_params = config['bakta_params']
    threads: 15
    log: analysis_dir + "/logs/pangenome/bakta_v3.log"
    conda: "bakta"
    message:
        "Running bakta to annotate the pangenome reference at ==>> {input}"
    shell:
        "bakta --db {params.bakta_params[db]} {params.bakta_params[gram]} --skip-trna --keep-contig-headers --genus Borrelia --species burgdorferi --strain pan-sensu-stricto --threads {threads} {params.bakta_params[opts]} --force --output {output.outdir} --prefix roary_pangene_v3 {input} 2>&1 >{log}"
