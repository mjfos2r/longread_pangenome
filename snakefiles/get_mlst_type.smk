import pandas
import glob
import mjf.tools as mjf

configfile: local("configs/singularity_config.yaml")
analysis_dir = config['analysis_dir']

rule all:
    input:
        analysis_dir+"/reports/mlst/mlst_results.allele.csv"

##
# MLST Type each assembly!
## 

rule getMLST_type:
    input: glob.glob(analysis_dir + "/assemblies/*/contigs/*.fasta",sample=SAMPLES)
    output:
        outdir = directory(analysis_dir+"/reports/mlst/"),
        final_out = analysis_dir+"/reports/mlst/mlst_results.allele.csv"
    singularity: "singularity/images/mlst_check_latest.sif"
    params:
        singularity_out =  analysis_dir+"/singularity_output/",
    threads: 360
    shell:
        "get_sequence_type -d {threads} -s 'Borrelia' -o {output.outdir} {input}"
