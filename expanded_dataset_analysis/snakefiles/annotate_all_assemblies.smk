import glob

configfile: "configs/pangenome_analysis.yaml"
analysis_dir = config['analysis_dir']

samples = {
    "shortread" : [s.split("/")[-1].split(".")[0] for s in glob.glob(analysis_dir+"/paired_assemblies/paired_only/shortread/contigs/*.fasta")],
    "longread" :  [s.split("/")[-1].split(".")[0] for s in glob.glob(analysis_dir+"/paired_assemblies/paired_only/longread/contigs/*.fasta")],
}

rule all:
    input:
        expand(analysis_dir+"/paired_assemblies/paired_only/{method}/annotation/{sample}/{sample}.gff3", method=["shortread"], sample=samples["shortread"]),
        expand(analysis_dir+"/paired_assemblies/paired_only/{method}/annotation/{sample}/{sample}.gff3", method=["longread"], sample=samples["longread"])

rule bakta:
    input:
        assembly = analysis_dir+"/paired_assemblies/paired_only/{method}/contigs/{sample}.fasta"
    output:
        outdir = directory(analysis_dir+"/paired_assemblies/paired_only/{method}/annotation/{sample}"),
        outgff = analysis_dir+"/paired_assemblies/paired_only/{method}/annotation/{sample}/{sample}.gff3",
        checkpoint = touch(analysis_dir + "/checkpoints/paired_only/{method}/.{sample}_bakta_finished")
    params:
        bakta_params = config['bakta_params']
    threads: 60
    log: analysis_dir + "/logs/{method}/{sample}.bakta.log"
    conda: "bakta"
    message:
        "Running bakta for sample {wildcards.sample}"
    shell:
        "bakta --db {params.bakta_params[db]} {params.bakta_params[gram]} --threads {threads} {params.bakta_params[opts]} --force --output {output.outdir} --prefix {wildcards.sample} {input.assembly} 2>&1 >{log}"

# cmd is:
# bakta --db /home/mf019/db/bakta --gram - --threads 60 --verbose --keep-contig-headers --force --output longread_analysis/paired_assemblies/paired_only/shortread/annotation/UNY172 --prefix UNY172 longread_analysis/paired_assemblies/paired_only/shortread/contigs/UNY172.fasta 2>&1 >longread_analysis/logs/shortread/UNY172.bakta.log
#rule Quast:
#    input:
#        input_assembly = analysis_dir + "/paired_assemblies/{method}/contigs/{sample}.fasta"
#    output:
#        outdir = directory(analysis_dir+"/reports/quast/{method}/{sample}/"),
#        report = analysis_dir+"/reports/quast/{method}/{sample}/report.txt",
#        checkpoint = touch(analysis_dir + "/checkpoints/{method}/.{sample}_quast_finished")
#    params:
#        quast_params = config['quast_params'],
#        ref_gff = config['reference_gff'],
#        ref_fa = config['reference_fa']
#    threads: 60
#    log: analysis_dir + "/logs/{method}/{sample}.quast.log"
#    conda: "quast"
#    message:
#        "Running Quast for {wildcards.method}, sample: {wildcards.sample}"
#    shell:
#        "quast -t {threads} {input.input_assembly} -r {params.ref_fa} -g {params.ref_gff} {params.quast_params[mode]} {params.quast_params[options]} "
#        "-o {output.outdir} 2>&1 >{log}"

#def all_quast_reps_exist(wildcards):
#    all_reps = []
#    for key,values in samples.items():
#        for value in values:
#            checkpoint_file = analysis_dir + f"/checkpoints/{key}/.{value}_quast_finished"
#            all_reps.append(checkpoint_file)
#    return all_reps

## Checkpoint rule to ensure all files exist
#checkpoint all_quast_reps_exist:
#    input:
#        all_quast_reps_exist
#    output:
#        touch(analysis_dir + "/checkpoints/.all_quast_reps_exist")
#

## Define input function to check if all sample BAM files exist
#def all_annotations_exist(wildcards):
#    all_reps = []
#    for key,values in samples.items():
#        for value in values:
#            checkpoint_file = analysis_dir + f"/checkpoints/{key}/.{value}_bakta_finished"
#            all_reps.append(checkpoint_file)
#    return all_reps
#
## Checkpoint rule to ensure all files exist
#checkpoint all_annotations_exist:
#    input:
#        all_annotations_exist
#    output:
#        touch(analysis_dir + "/checkpoints/.all_annotations_exist")
#

#rule agat_basic_stats:
#    input:
#        analysis_dir+"assemblies/{method}/annotation/{sample}.gff3"
#    output:
#        analysis_dir+"/reports/annotation/{method}/{sample}.basic_stats.txt"
#    conda: "bakta"
#    log: analysis_dir+"/reports/annotation/{method}/{sample}.basic_stats.log"
#    message:
#        "Running agat_basic_stats for sample {wildcards.sample}"
#    shell:
#        "agat_sq_stat_basic.pl --gff {input} --output {output} 2>{log}"
#
#rule agat_long_stats:
#    input:
#        analysis_dir+"assemblies/{method}/annotation/{sample}.gff3"
#    output:
#        analysis_dir+"/reports/annotation/{method}/{sample}.long_stats.txt"
#    conda: "bakta"
#    log: analysis_dir+"/reports/annotation/{method}/{sample}.long_stats.log"
#    message:
#        "Running agat_long_stats for sample {wildcards.sample}"
#    shell:
#        "agat_sp_statistics.pl --gff {input} --output {output} 2>{log}"
#
#rule MultiQC:
#    input:
#        checkpoint_quast = analysis_dir + "/checkpoints/.all_quast_reps_exist",
#        checkpoint_bakta = analysis_dir + "/checkpoints/.all_annotations_exist",
#    output:
#        outdir = directory(analysis_dir + "/reports/multiQC/"),
#        outfiles = analysis_dir + "/reports/multiQC/multiqc_report.html"
#    conda: "multiqc"
#    threads: 360
#    params:
#        reports = analysis_dir,
#        ignore = "--ignore '*.conf' --ignore '*/tmp/*'"
#    message: "Running MultiQC to compile all reports into a single html document!"
#    shell:
#        "multiqc --interactive {params.reports} --outdir {output.outdir} {params.reports} {params.ignore}"
