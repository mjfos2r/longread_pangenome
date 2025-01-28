import pandas
import glob

configfile: local("configs/pangenome_analysis.yaml")
analysis_dir = config['analysis_dir']
version = config['version']

<<<<<<< HEAD
ASSEMBLIES = [os.path.basename(os.path.dirname(f)) for f in glob.glob(f"{analysis_dir}/assemblies/*/*.fna")]
print(ASSEMBLIES)
rule all:
    input:
        expand(f'{analysis_dir}/kraken2/reports/{{assembly}}_kraken2.report', assembly=ASSEMBLIES),
        expand(f'{analysis_dir}/kraken2/reports/{{assembly}}_kraken2.output', assembly=ASSEMBLIES),
=======
ASSEMBLIES = [os.path.basename(os.path.dirname(f)) for f in glob.glob(f"{analysis_dir}/assemblies/expanded_v2/*/*.fna")]
print(ASSEMBLIES)
rule all:
    input:
        expand(f'{analysis_dir}/kraken2/reports/{{assembly}}.report', assembly=ASSEMBLIES),
        expand(f'{analysis_dir}/kraken2/reports/{{assembly}}.output', assembly=ASSEMBLIES),
>>>>>>> ce4e92d (this will break everything)
        f'{analysis_dir}/kraken2/kraken2_report.krona',
        f'{analysis_dir}/kraken2/kraken2_report.krona.html'

rule run_Kraken2:
<<<<<<< HEAD
    input: f'{analysis_dir}/assemblies/{{assembly}}/{{assembly}}.fna'
    output:
        output = f'{analysis_dir}/kraken2/reports/{{assembly}}_kraken2.output',
        report = f'{analysis_dir}/kraken2/reports/{{assembly}}_kraken2.report'
        #outdir = directory(f'{analysis_dir}/kraken2/reports/{{assembly}}/')
    singularity: 'singularity/images/kraken2_latest.sif'
    log: f'{analysis_dir}/logs/kraken2/kraken2.log'
=======
    input: f'{analysis_dir}/assemblies/expanded_v2/{{assembly}}/{{assembly}}.fna'
    output:
        output = f'{analysis_dir}/kraken2/reports/{{assembly}}.output',
        report = f'{analysis_dir}/kraken2/reports/{{assembly}}.report'
        #outdir = directory(f'{analysis_dir}/kraken2/reports/{{assembly}}/')
    singularity: 'singularity/images/kraken2_latest.sif'
    log: f'{analysis_dir}/logs/kraken2/{{assembly}}_kraken2.log'
>>>>>>> ce4e92d (this will break everything)
    params:
        db = config['kraken2_db']
    threads: 32
    shell:
        """
<<<<<<< HEAD
        (kraken2 --use-names --memory-mapping --db {params.db} \
            --threads {threads} \
            --output {output.output} \
            --report {output.report} \
            {input}) 2&1 >{log}
=======
        kraken2 --use-names --memory-mapping --db /kraken2_dbs \
            --threads {threads} \
            --output {output.output} \
            --report {output.report} \
            {input}
>>>>>>> ce4e92d (this will break everything)
        """

rule krona_plot:
    input:
<<<<<<< HEAD
        report = f'{analysis_dir}/kraken2/reports/{{assembly}}_kraken2.report'
    output:
        krona = f'{analysis_dir}/kraken2/krona/{{assembly}}_kraken2_report.krona',
    singularity: 'singularity/images/kraken2_latest.sif'
=======
        report = f'{analysis_dir}/kraken2/reports/{{assembly}}.report'
    output:
        krona = f'{analysis_dir}/kraken2/krona/{{assembly}}.krona',
    singularity: 'singularity/images/kraken2_latest.sif'
    threads: 8
>>>>>>> ce4e92d (this will break everything)
    shell:
        """
        kreport2krona.py -r {input} -o {output.krona}
        """

rule html_plot:
    input:
<<<<<<< HEAD
        krona = expand(f'{analysis_dir}/kraken2/krona/{{assembly}}_kraken2_report.krona', assembly=ASSEMBLIES)
    output:
        html = f'{analysis_dir}/kraken2/kraken2_report.krona.html'
    singularity: 'singularity/images/kraken2_latest.sif'
=======
        krona = expand(f'{analysis_dir}/kraken2/krona/{{assembly}}.krona', assembly=ASSEMBLIES)
    output:
        html = f'{analysis_dir}/kraken2/kraken2_report.krona.html',
        report = f'{analysis_dir}/kraken2/kraken2_report.krona'
    singularity: 'singularity/images/kraken2_latest.sif'
    threads: 8
>>>>>>> ce4e92d (this will break everything)
    params:
        r_files = lambda wildcards, input: " ".join([""] + input.krona)
    shell:
        """
        ktImportText {params.r_files} -o {output.html}
        """