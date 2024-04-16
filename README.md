# Working Repository for our Borrelia burgdorferi Genome Wide Association Study
This is the repository for the longread GWAS project performed by the Lemieux Lab at MGH in Boston.
Snakefiles are in various states of functionality (issues stem from directory clobbering and this is a known limitation of snakemake workflows)
However, core functionality is outlined in each rule should the snakefile fail to successfully complete.

`longread_GWAS/.git`
- git repo stuff goes here, don't touch unless you absolutely need to.

`longread_GWAS/.github`
- github actions and automations go here (leftover from the snakemake workflow template I originally forked for this repo, should probably tidy up)

`longread_GWAS/.vscode`
- metadata and preferences for the vscode project are stored here.

`longread_GWAS/assemblies`
 - here is where the assembly fastas and annotations in GFF3 & GBFF are housed.
    within this, there are subdirectories for hybrid assemblies, illumina assemblies, nanopore assemblies, and pacbio assemblies. Isolates with paired assemblies across multiple methods are housed in the `longread_GWAS/assemblies/paired_assemblies` directory.
    it is further subdivided into `annotation/`, `contigs/`, `shortread/` and `paired_only/`. the first two are for all paired assemblies together which includes longread assemblies without a shortread pair. `shortread/` is just illumina assemblies from this subset, and `paired_only/` is broken down between short and long read assemblies and further between `annotation/` and `contigs/`.
    {{TODO}} deal with massive redundancy for storage optimization.

`longread_GWAS/configs`
- configs for the snakemake workflows are kept here. see snakefiles below for further details.

`longread_GWAS/longread_analysis`
- this is the working and output directory for the snakemake workflows defined in snakefiles.
    due to file size, this is UNTRACKED by github and this needs to have bucket backup automated at some point.

`longread_GWAS/notebooks`
- jupyter notebooks primarily for metadata wranglin' are kept here.
- jupyter notebooks and data for OspC and MLST typing are also kept here.

`longread_GWAS/plasmid_id`
- this is where notebooks and associated data are kept for plasmid abundance analysis.

`longread_GWAS/ref`
- reference genome(s) and annotations are kept here.

`longread_GWAS/scripts`
- assorted scripts are housed here. This is where scripts for use in snakemake workflows are kept.
{{TODO}} this needs to be cleaned up and further populated.

`longread_GWAS/singularity`
- here is where container defs/recipes and container images are stored. This project uses Singularity/apptainer containers *EXCLUSIVELY*.
any docker container will necessarily be converted into an apptainer for integration with snakemake.

`longread_GWAS/snakefiles`
- here is where the main snakefiles for the workflows are housed.
configs for these workflows are kept in `longread_GWAS/configs` and these configs define paths, working directories, and rule params.

`longread_GWAS/spreadsheets`
- here is where unorganized/assorted spreadsheets are housed. {{TODO}} these need to be sorted into their corresponding subdirectories elsewhere.

{{TODO}}
1. organize and optimize this project's repo to streamline further analyses and development
2. add in submodule functionality via git to incorporate custom python packages I've written elsewhere.
    - lipoPredict
    - plasmidBlaster
    - .....?
3. move everything to singularity containers for scalability
4. K8s implementation?
5. Fix snakemake workflows involving Roary since pathnames are clobbered by snakemake's output directory creation.
    - for `rule Roary:`, `output: dir(path/to/output/{roary_output})` is required to satisfy DAG creation in the rest of the workflow, however, this involves snakemake creating output directories prior to rule execution. Roary, when given an output directory via CLI flags, does not like overwriting directories so appends a timestamp at the end to differentiate different roary runs using the same run_ID. This prevents downstream rules from locating the roary output files and crashes the workflow.
    This could likely be dealt with using a checkpoint function but I honestly cannot bring myself to write one. May just work around this and write a simple rule to `mv` all of the roary output to the directory that snakemake expects. // alternatively, I could just write a simple bash script and call that with snakemake? though I don't think that solves my issue of snakemake creating the output directory before rule execution c: ....