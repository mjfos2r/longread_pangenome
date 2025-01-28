import glob
import os

configfile: "configs/pangenome_analysis.yaml"
analysis_dir = config['analysis_dir']
input_dir = config['input_dir']

# Assuming directory structure like:
# input_dir/
#   genome1/
#     genome1.fna
#   genome2/
#     genome2.fnaf

GENOMES = [os.path.basename(os.path.dirname(f)) for f in glob.glob(f"{input_dir}/*/*.fna")]

<<<<<<< HEAD
rule all:
    input:
        [f"{analysis_dir}/assemblies/{genome}/{genome}.gff3" for genome in GENOMES],
        [f"{analysis_dir}/bakta/{genome}.long_stats.txt" for genome in GENOMES],
        [f"{analysis_dir}/bakta/{genome}.basic_stats.txt" for genome in GENOMES],
        [f"{analysis_dir}/quast/{genome}/report.txt" for genome in GENOMES],
        f"{analysis_dir}/multiQC/multiqc_report.html"

rule bakta:
    input:
        assembly = f"{input_dir}/{{genome}}/{{genome}}.fna",
        replicons = f"{input_dir}/{{genome}}/{{genome}}.contig_names.tsv"
    output:
        outdir = directory(f"{analysis_dir}/assemblies/{{genome}}/"),
        outgff = f"{analysis_dir}/assemblies/{{genome}}/{{genome}}.gff3"
    params:
        bakta_params = config['bakta_params']
    threads: 32
    log: f"{analysis_dir}/logs/bakta/{{genome}}.bakta.log"
    singularity: "singularity/images/bakta_1.8.2--pyhdfd78af_0.sif"
    message:
        "Running bakta for sample {wildcards.genome}"
    shell:
        "bakta --db {params.bakta_params[db]} {params.bakta_params[gram]} --threads {threads} "
        "{params.bakta_params[opts]} --force --output {output.outdir} "
        "--prefix {wildcards.genome} --replicons {input.replicons} {input.assembly} 2>&1 >{log}"
=======

# Debug prints for expected files
print("\nExpected files for rule all:")
print("\nGFF files:")
print([f"{analysis_dir}/assemblies/{genome}/{genome}.gff3" for genome in GENOMES])
print("\nFNA files:")
print([f"{analysis_dir}/assemblies/{genome}/{genome}.fna" for genome in GENOMES])
print("\nLong stats:")
print([f"{analysis_dir}/reports/agat/{genome}.long_stats.txt" for genome in GENOMES])
print("\nBasic stats:")
print([f"{analysis_dir}/reports/agat/{genome}.basic_stats.txt" for genome in GENOMES])
print("\nQuast reports:")
print([f"{analysis_dir}/reports/quast/{genome}/report.txt" for genome in GENOMES])

# This is required because the stupid snakes cannot resolve these wildcards without having a major stroke and falling into
# an infinite recursion during DAG resolution. I cannot put into words how much I absolutely detest snakemake.
wildcard_constraints:
    genome="[^/]+"

rule all:
    input:
        [f"{analysis_dir}/assemblies/{genome}/{genome}.gff3" for genome in GENOMES],
        [f"{analysis_dir}/assemblies/{genome}/{genome}.fna" for genome in GENOMES],
        [f"{analysis_dir}/reports/agat/{genome}.long_stats.txt" for genome in GENOMES],
        [f"{analysis_dir}/reports/agat/{genome}.basic_stats.txt" for genome in GENOMES],
        [f"{analysis_dir}/reports/quast/{genome}/report.txt" for genome in GENOMES],
        f"{analysis_dir}/reports/multiQC/multiqc_report.html"

rule bakta:
    input:
        assembly = f"{input_dir}/{{genome}}/{{genome}}.fna"
    output:
        outdir = directory(f"{analysis_dir}/assemblies/{{genome}}/"),
        outgff = f"{analysis_dir}/assemblies/{{genome}}/{{genome}}.gff3",
        outfna = f"{analysis_dir}/assemblies/{{genome}}/{{genome}}.fna"
    params:
        gram = config['bakta']['gram'],
        opts = config['bakta']['opts']
    threads: 64
    log: f"{analysis_dir}/logs/bakta/{{genome}}.bakta.log"
    singularity: "singularity/images/bakta_latest.sif"
    message:
        "Running bakta for sample {wildcards.genome}"
    shell:
        "(bakta --db /db {params.gram} --threads {threads} "
        "{params.opts} --output {output.outdir} "
        "--prefix {wildcards.genome} --keep-contig-headers -v --skip-plot {input.assembly}) 2>&1 >{log} &&"
        "sed -i -E '/^##feature-ontology/s|http[s]?://[^ ]+|https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/refs/tags/v3.1/so.obo|' {output.outgff}"
>>>>>>> ce4e92d (this will break everything)

rule Quast:
    input:
        f"{analysis_dir}/assemblies/{{genome}}/{{genome}}.fna"
    output:
<<<<<<< HEAD
        outdir = directory(f"{analysis_dir}/quast/{{genome}}"),
        report = f"{analysis_dir}/quast/{{genome}}/report.txt"
    params:
        quast_params = config['quast_params'],
        ref_gff = config['reference_gff'],
        ref_fa = config['reference_fa']
    threads: 32
=======
        outdir = directory(f"{analysis_dir}/reports/quast/{{genome}}"),
        report = f"{analysis_dir}/reports/quast/{{genome}}/report.txt"
    params:
        mode = config['quast']['mode'],
        opts = config['quast']['opts'],
        ref_gff = config['reference_gff'],
        ref_fa = config['reference_fa']
    threads: 64
>>>>>>> ce4e92d (this will break everything)
    log: f"{analysis_dir}/logs/quast/{{genome}}.quast.log"
    singularity: "singularity/images/quast_latest.sif"
    message:
        "Running Quast for {wildcards.genome}"
    shell:
<<<<<<< HEAD
        "quast -t {threads} {input} -r {params.ref_fa} -g {params.ref_gff} {params.quast_params[mode]} "
        "{params.quast_params[options]} -o {output.outdir} 2>&1 >{log}"
=======
        "quast.py -t {threads} {input} -r {params.ref_fa} -g {params.ref_gff} {params.mode} "
        "{params.opts} -o {output.outdir} 2>&1 >{log}"
>>>>>>> ce4e92d (this will break everything)

rule agat_basic_stats:
    input:
        f"{analysis_dir}/assemblies/{{genome}}/{{genome}}.gff3"
    output:
<<<<<<< HEAD
        f"{analysis_dir}/bakta/{{genome}}.basic_stats.txt"
    singularity: "/home/mf019/longread_pangenome/singularity/images/agat_1.4.1--pl5321hdfd78af_0.sif"
    threads: 16
    log: f"{analysis_dir}/logs/bakta/{{genome}}.basic_agat.log"
=======
        f"{analysis_dir}/reports/agat/{{genome}}.basic_stats.txt"
    singularity: "/home/mf019/longread_pangenome/singularity/images/agat_fixed.sif"
    threads: 64
    retries: 3
    log: f"{analysis_dir}/logs/agat/{{genome}}.basic_agat.log"
>>>>>>> ce4e92d (this will break everything)
    message:
        "Running agat_basic_stats for assembly: {wildcards.genome}"
    shell:
        "agat_sq_stat_basic.pl --gff {input} --output {output} 2>&1 >{log}"

rule agat_long_stats:
    input:
        f"{analysis_dir}/assemblies/{{genome}}/{{genome}}.gff3"
    output:
<<<<<<< HEAD
        f"{analysis_dir}/bakta/{{genome}}.long_stats.txt"
    singularity: "/home/mf019/longread_pangenome/singularity/images/agat_1.4.1--pl5321hdfd78af_0.sif"
    threads: 16
    log: f"{analysis_dir}/logs/bakta/{{genome}}.long_agat.log"
=======
        f"{analysis_dir}/reports/agat/{{genome}}.long_stats.txt"
    singularity: "/home/mf019/longread_pangenome/singularity/images/agat_fixed.sif"
    threads: 64
    retries: 3
    log: f"{analysis_dir}/logs/agat/{{genome}}.long_agat.log"
>>>>>>> ce4e92d (this will break everything)
    message:
        "Running agat_long_stats for assembly: {wildcards.genome}"
    shell:
        "agat_sp_statistics.pl --gff {input} --output {output} 2>&1 >{log}"

rule MultiQC:
<<<<<<< HEAD
    input:
        rules.all.input[:-1]  # All inputs from rule all except the MultiQC report itself
    output:
        outdir = directory(f"{analysis_dir}/annotation/multiQC/"),
        outfiles = f"{analysis_dir}/multiQC/multiqc_report.html"
    singularity: "singularity/images/multiqc_latest.sif"
    threads: 32
    params:
        reports = analysis_dir,
        ignore = "--ignore '*.conf' --ignore '*/tmp/* --ignore ncbi_dataset/*'"
    message: "Running MultiQC to compile all reports into a single html document!"
    shell:
        "multiqc --interactive --outdir {output.outdir} {params.reports} {params.ignore}"
=======
   input:
       rules.all.input[:-1]
   output:
       outdir = directory(f"{analysis_dir}/reports/multiQC/"),
       outfiles = f"{analysis_dir}/reports/multiQC/multiqc_report.html"
   singularity: "singularity/images/multiqc_latest.sif"
   threads: 64
   params:
       input_dir = analysis_dir,
       ignore = "--ignore '*.conf' --ignore '*/tmp/*' --ignore 'ncbi_dataset/*'"
   message: "Running MultiQC to compile all reports into a single html document!"
   shell:
       """
       multiqc --interactive \
       --outdir {output.outdir} \
       {params.input_dir} \
       {params.ignore} \
       -v
       """
>>>>>>> ce4e92d (this will break everything)
