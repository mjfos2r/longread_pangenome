import pandas
import glob
import mjf.tools as mjf

configfile: local("configs/singularity_config.yaml")
analysis_dir = config['analysis_dir']

rule all:
    input:
        analysis_dir+"/reports/mlst/mlst_results.allele.csv"

samples = {"illumina" : [v.split("/")[-1].split(".")[0] for v in glob.glob(analysis_dir + "/assemblies/illumina/contigs/*.fasta")],
             "pacbio" : [s.split("/")[-1].split(".")[0] for s in glob.glob(analysis_dir + "/assemblies/pacbio/contigs/*.fasta")],
           "nanopore" : [t.split("/")[-1].split(".")[0] for t in glob.glob(analysis_dir + "/assemblies/nanopore/contigs/*.fasta")],
             "hybrid" : [u.split("/")[-1].replace(".fasta", "") for u in glob.glob(analysis_dir + "/assemblies/hybrid/contigs/*.fasta")]}

all_input = []
for key,values in samples.items():
    for value in values:
        all_input.append(f"{analysis_dir}/assemblies/{method}/contigs/{sample}.fasta")

##
# MLST Type each assembly!
## 

rule getMLST_type:
    input: expand(f"{analysis_dir}/assemblies/{method}/contigs/{sample}.fasta", analysis_dir=analysis_dir, method=samples.keys(), sample=samples.items())
    output:
        outdir = directory(analysis_dir+"/reports/mlst/"),
        final_out = analysis_dir+"/reports/mlst/mlst_results.allele.csv"
    singularity: "singularity/images/mlst_check_latest.sif"
    params:
        singularity_out =  analysis_dir+"/singularity_output/",
    threads: 360
    message: "{input}"
    shell:
        "get_sequence_type -d {threads} -s 'Borrelia' -o {output.outdir} {input}"
