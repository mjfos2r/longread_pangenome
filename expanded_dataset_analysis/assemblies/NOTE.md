# Notes on redundant GFF directory

Due to how Roary runs within snakemake, for whatever reason, I have not been able to get it to correctly accept inputs from a glob. 
For whatever reason, it works best to just copy all of the gff3 files to be analyzed into their own directory and feed that as input.
This is not optimal.

The code to perform this is trivial.
inside of this directory (assemblies), run the following command (use whatever version you want if you're running on new data):
```{bash}
mkdir gffs_v5 && cp dataset_v5/*/*.gff3 gffs_v5/
```

To check to make sure everything is groovy:

```{bash}
echo "All GFFs: $( ls -1 dataset_v5/*/*.gff3 | wc -l )"
echo "Copied GFFs: $( ls -1 gffs_v5/*.gff3 | wc -l )"
```