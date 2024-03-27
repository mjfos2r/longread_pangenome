import glob

configfile: 'configs/pangenome_analysis.yaml'

analysis_dir = config['analysis_dir']

rule all:
    input:
        #f'{analysis_dir}/pangenome/v3/roary_longread_v1/core_gene_alignment.aln',
        f'{analysis_dir}/pangenome/v3/tree/longread_v1/roary_longread_v1_core_aln.newick',
        f'{analysis_dir}/pangenome/v3/tree/longread_v1/RAxML_bestTree.roary_longread_v1_gtrgam',
        f'{analysis_dir}/pangenome/v3/bakta_longread_v1/roary_core_longread_v1.gff3'

#rule Roary:
#    input:
#        glob.glob(f'{analysis_dir}/paired_assemblies/annotation/longread/*.gff3')
#    output:
#        outdir = directory(f'{analysis_dir}/pangenome/v3/roary_longread_v1/'),
#        outfile = f'{analysis_dir}/pangenome/v3/roary_longread_v1/core_gene_alignment.aln',
#        outasm = f'{analysis_dir}/pangenome/v3/roary_longread_v1/pan_genome_reference.fa',
#    threads: 15
#    message: "running roary on all of our annotations! : \n roary -p {threads} -f {output.outdir} {input} 2>&1 >{log}"
#    log: f'{analysis_dir}/pangenome/v3/logs/roary_log.txt'
#    singularity: 'singularity/images/roary_latest.sif'
#    shell:
#        "roary -e --mafft -p {threads} -f {output.outdir} {input} >{log} 2>&1"

rule FastTree:
    input:
        f'{analysis_dir}/pangenome/v3/roary_longread_v1/core_gene_alignment.aln'
    output:
        f'{analysis_dir}/pangenome/v3/tree/longread_v1/roary_longread_v1_core_aln.newick'
    threads: 15
    conda: "fasttree"
    message: "running fasttree on the roary output ==>> {input}"
    log: f'{analysis_dir}/pangenome/v3/logs/fasttree.txt'
    shell:
        "FastTreeMP -gamma -nt -gtr {input} > {output} 2>{log}"

rule RAxML:
    input:
        f'{analysis_dir}/pangenome/v3/roary_longread_v1/core_gene_alignment.aln'
    output:
        outdir = directory(f'{analysis_dir}/pangenome/v3/tree/longread_v1/'),
        outfile = f'{analysis_dir}/pangenome/v3/tree/longread_v1/RAxML_bestTree.roary_longread_v1_gtrgam'
    threads: 15
    conda: "raxml"
    message: "running RAxML on the roary output ==>> {input}"
    log: f'{analysis_dir}/pangenome/v3/logs/RAxML_log.txt'
    shell:
        "raxmlHPC -m GTRGAMMA -p 71235 -s {input} -n roary_longread_v1_gtrgam -w {output.outdir} >{log} 2>&1"

rule bakta:
    input: f'{analysis_dir}/pangenome/v3/roary_longread_v1/pan_genome_reference.fa'
    output:
        outdir = directory(f'{analysis_dir}/pangenome/v3/bakta_longread_v1/'),
        outgff = f'{analysis_dir}/pangenome/v3/bakta_longread_v1/roary_core_longread_v1.gff3',
    params:
        bakta_params = config['bakta_params'],
        bakta_prefix = 'longread_v1'

    threads: 15
    log: f'{analysis_dir}/pangenome/v3/logs/bakta_longread_v1.log'
    conda: "bakta"
    message:
        "Running bakta to annotate the pangenome reference at ==>> {input}"
    shell:
        "bakta --db={params.bakta_params[db]} {params.bakta_params[gram]} --skip-trna --keep-contig-headers --genus Borrelia --species burgdorferi --strain roary_pangenome_v3 --threads {threads} {params.bakta_params[opts]} --force --output {output.outdir} --prefix {params.bakta_prefix} {input} 2>&1 >{log}"