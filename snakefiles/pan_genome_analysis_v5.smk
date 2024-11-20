import glob

configfile: 'configs/pangenome_analysis.yaml'

analysis_dir = config['analysis_dir']
gffs_dir = config['gffs_dir']
version = config['version']

rule all:
    input:
        f'{analysis_dir}/results/{version}/roary_{version}/core_gene_alignment.aln',
        f'{analysis_dir}/results/{version}/tree/{version}/fasttree_roary_core_v5.newick',
        f'{analysis_dir}/results/{version}/tree/{version}/RAxML_bestTree.roary_v5_1_gtrgam',
        f'{analysis_dir}/results/{version}/bakta_{version}/roary_core_v5.gff3'

rule Roary:
    input:
        glob.glob(f'{gffs_dir}/*.gff3')
    output:
        outfile = f'{analysis_dir}/results/{version}/roary_{version}/core_gene_alignment.aln',
        outasm = f'{analysis_dir}/results/{version}/roary_{version}/pan_genome_reference.fa',
    params:
        anal = f'{analysis_dir}/results/v5',
        outdir = f'{analysis_dir}/results/{version}/roary_v5',
        temp_dir = f'{analysis_dir}/results/{version}/roary_v5_temp'
    threads: 128
    message: "running roary on all of our annotations! : \n roary -p {threads} -f {params.outdir} {input} 2>&1 >{log}"
    log: f'{analysis_dir}/results/{version}/logs/roary_log.txt'
    singularity: 'singularity/images/roary_latest.sif'
    shell: # -s to not split paralogs
        """
        roary -s -e --mafft -p {threads} -f {params.temp_dir} {input} >{log} 2>&1
        # move results from the roary_temp dir to the snakemake dir.
        # I hate snakes
        ls {params.anal}
        mv {params.temp_dir}/* {params.outdir}/

        """

rule FastTree:
    input:
        f'{analysis_dir}/results/{version}/roary_{version}/core_gene_alignment.aln'
    output:
        f'{analysis_dir}/results/{version}/tree/{version}/fasttree_roary_core_v5.newick'
    threads: 128
    singularity: 'singularity/images/fasttree_latest.sif'
    message: "running fasttree on the roary output ==>> {input}"
    log: f'{analysis_dir}/results/{version}/logs/fasttree.txt'
    shell:
        "FastTree -gamma -nt -gtr {input} > {output} 2>{log}"

rule RAxML:
    input:
        f'{analysis_dir}/results/{version}/roary_{version}/core_gene_alignment.aln'
    output:
        outdir = directory(f'{analysis_dir}/results/{version}/tree/{version}/'),
        outfile = f'{analysis_dir}/results/{version}/tree/{version}/RAxML_bestTree.roary_v5_1_gtrgam'
    threads: 128
    singularity: 'singularity/images/raxml_latest.sif'
    message: "running RAxML on the roary output ==>> {input}"
    log: f'{analysis_dir}/results/{version}/logs/RAxML_log.txt'
    shell:
        "raxmlHPC -m GTRGAMMA -p 71235 -s {input} -n roary_v5_gtrgam -w {output.outdir} >{log} 2>&1"

rule bakta:
    input: f'{analysis_dir}/results/{version}/roary_{version}/pan_genome_reference.fa'
    output:
        outdir = directory(f'{analysis_dir}/results/{version}/bakta_{version}/'),
        outgff = f'{analysis_dir}/results/{version}/bakta_{version}/roary_core_v5.gff3',
    params:
        gram = config['bakta']['gram'],
        opts = config['bakta']['opts']
    threads: 128
    log: f'{analysis_dir}/results/{version}/logs/bakta_v5.log'
    singularity: "singularity/images/bakta_latest.sif"
    message:
        "Running bakta for pangenome!"
    shell:
        "(bakta --db /db {params.gram} --threads {threads} "
        "{params.opts} --output {output.outdir} "
        "--prefix roary_pangenome_v5 --genus Borrelia --species burgdorferi --strain sensu-stricto --keep-contig-headers -v --skip-plot {input}) 2>&1 >{log} &&"
        "sed -i -E '/^##feature-ontology/s|http[s]?://[^ ]+|https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/refs/tags/v3.1/so.obo|' {output.outgff}"
