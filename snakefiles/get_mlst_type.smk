# import pandas to parse the sample file
import pandas
import mjf.tools as mjf

# Load the config file
configfile: local("configs/singularity_config.yaml")

analysis_dir = config['analysis_dir']
sample_sheet = config['sample_sheet']

SAMPLES = mjf.parseSampleSheet(sample_sheet)
assemblers = ['flye','canu']

rule all:
    input:
        expand(analysis_dir+"/reports/assembly/denovo-telomere/{assembler}/mlst_results.allele.csv",assembler=assemblers)

##
# MLST Type each assembly!
## 

rule getMLST_type:
    input: expand(analysis_dir+"/samples/{sample}/assembly/telomeres/{{assembler}}/polished/{sample}.consensus.fasta",sample=SAMPLES),
    output: analysis_dir+"/reports/assembly/denovo-telomere/{assembler}/mlst_results.allele.csv"
    conda: mlst
    singularity: "singularity/mlst_check_latest.sif"
    threads: 360
    params: params
    message: "getting MLST type for {wildcards.assembler} generated assemblies!"
    shell:
        "export MLST_DATABASES=/home/db/mlst/ &&"
        "get_sequence_type -s 'Borrelia' {input}"
