import glob

configfile: "configs/annotate_all_assemblies.yaml"
analysis_dir = config['analysis_dir']

pacbio_samples = [s.split("/")[-1].split(".")[0] for s in glob.glob(analysis_dir + "/assemblies/pacbio/contigs/*.fasta")]
nanopore_samples = [t.split("/")[-1].split(".")[0] for t in glob.glob(analysis_dir + "/assemblies/nanopore/contigs/*.fasta")]
hybrid_samples = [u.split("/")[-1].replace(".fasta", "") for u in glob.glob(analysis_dir + "/assemblies/hybrid/contigs/*.fasta")]
illumina_samples = [v.split("/")[-1].split(".")[0] for v in glob.glob(analysis_dir + "/assemblies/illumina/contigs/*.fasta")]

rule all:
    input:
        analysis_dir + "/reports/multiQC/multiqc_report.html",
        expand(analysis_dir+"/assemblies/pacbio/annotation/{pb_sample}/{pb_sample}.gff3",pb_sample=pacbio_samples),
        expand(analysis_dir+"/reports/quast/pacbio/{pb_sample}/report.txt",pb_sample=pacbio_samples),
        expand(analysis_dir+"/assemblies/nanopore/annotation/{ont_sample}/{ont_sample}.gff3",ont_sample=nanopore_samples),
        expand(analysis_dir+"/reports/quast/nanopore/{ont_sample}/report.txt",ont_sample=nanopore_samples),
        expand(analysis_dir+"/assemblies/hybrid/annotation/{hy_sample}/{hy_sample}.gff3",hy_sample=hybrid_samples),
        expand(analysis_dir+"/reports/quast/hybrid/{hy_sample}/report.txt",hy_sample=hybrid_samples),
        expand(analysis_dir+"/assemblies/illumina/annotation/{il_sample}/{il_sample}.gff3",il_sample=illumina_samples),
        expand(analysis_dir+"/reports/quast/illumina/{il_sample}/report.txt",il_sample=illumina_samples)

rule Quast:
    input:
        analysis_dir + "/assemblies/{method}/contigs/{sample}.fasta",
    output:
        outdir = directory(analysis_dir+"/reports/quast/{method}/{sample}"),
        report = analysis_dir+"/reports/quast/{method}/{sample}/report.txt",
        checkpoint = touch(analysis_dir + "/checkpoints/{method}/.{sample}_quast_finished")
    params:
        quast_params = config['quast_params'],
        ref_gff = config['reference_gff'],
        ref_fa = config['reference_fa']
    threads: 60
    log: analysis_dir + "/logs/{method}/{sample}.quast.log"
    conda: "quast"
    message:
        "Running Quast for all samples"
    shell:
        "quast -t {threads} {input} -r {params.ref_fa} -g {params.ref_gff} {params.quast_params[mode]} {params.quast_params[options]} "
        "-o {output.outdir} 2>&1 >{log}"

def all_quast_reps_exist(wildcards):
    method_files = glob.glob(analysis_dir + "/assemblies/*/contigs/*.fasta")
    methods = [m.split("/")[-3] for m in method_files]
    samples = [m.split("/")[-1].replace(".fasta", "") for m in method_files]
    return expand(analysis_dir + "/checkpoints/{method}/.{sample}_quast_finished", method=methods, sample=samples)


# Checkpoint rule to ensure all files exist
checkpoint all_quast_reps_exist:
    input:
        all_quast_reps_exist
    output:
        touch(analysis_dir + "/checkpoints/.all_quast_reps_exist")

##
# annotate with bakta
##

rule bakta:
    input:
        assembly = analysis_dir+"/assemblies/{method}/contigs/{sample}.fasta"
    output:
        outdir = directory(analysis_dir+"/assemblies/{method}/annotation/{sample}"),
        outgff = analysis_dir+"/assemblies/{method}/annotation/{sample}/{sample}.gff3",
        checkpoint = touch(analysis_dir + "/checkpoints/{method}/.{sample}_bakta_finished")
    params:
        bakta_params = config['bakta_params']
    threads: 60
    log: analysis_dir + "/logs/{method}/{sample}.bakta.log"
    conda: "bakta"
    message:
        "Running bakta for sample {wildcards.sample}"
    shell:
        "bakta --db {params.bakta_params[db]} {params.bakta_params[gram]} --threads {threads} {params.bakta_params[opts]} --force --output {output.outdir} --prefix {wildcards.sample} {input.assembly} 2>&1 >{log}"

# Define input function to check if all sample BAM files exist
def all_annotations_exist(wildcards):
    method_files = glob.glob(analysis_dir + "/assemblies/*/contigs/*.fasta")
    methods = [m.split("/")[-3] for m in method_files]
    samples = [m.split("/")[-1].replace(".fasta", "") for m in method_files]
    return expand(analysis_dir + "/checkpoints/{method}/.{sample}_bakta_finished", method=methods, sample=samples)

# Checkpoint rule to ensure all files exist
checkpoint all_annotations_exist:
    input:
        all_annotations_exist
    output:
        touch(analysis_dir + "/checkpoints/.all_annotations_exist")


rule agat_basic_stats:
    input:
        analysis_dir+"assemblies/{method}/annotation/{sample}.gff3"
    output:
        analysis_dir+"/reports/annotation/{method}/{sample}.basic_stats.txt"
    conda: "bakta"
    log: analysis_dir+"/reports/annotation/{method}/{sample}.basic_stats.log"
    message:
        "Running agat_basic_stats for sample {wildcards.sample}"
    shell:
        "agat_sq_stat_basic.pl --gff {input} --output {output} 2>{log}"

rule agat_long_stats:
    input:
        analysis_dir+"assemblies/{method}/annotation/{sample}.gff3"
    output:
        analysis_dir+"/reports/annotation/{method}/{sample}.long_stats.txt"
    conda: "bakta"
    log: analysis_dir+"/reports/annotation/{method}/{sample}.long_stats.log"
    message:
        "Running agat_long_stats for sample {wildcards.sample}"
    shell:
        "agat_sp_statistics.pl --gff {input} --output {output} 2>{log}"

rule MultiQC:
    input:
        checkpoint_quast = analysis_dir + "/checkpoints/.all_quast_reps_exist",
        checkpoint_bakta = analysis_dir + "/checkpoints/.all_annotations_exist",
    output:
        outdir = directory(analysis_dir + "/reports/multiQC/"),
        outfiles = analysis_dir + "/reports/multiQC/multiqc_report.html"
    conda: "multiqc"
    threads: 360
    params: 
        reports = analysis_dir,
        ignore = "--ignore '*.conf' --ignore '*/tmp/*'"
    message: "Running MultiQC to compile all reports into a single html document!"
    shell:
        "multiqc --interactive {params.reports} --outdir {output.outdir} {params.reports} {params.ignore}"
