#!/usr/bin/env python3
#################################################
import os
import glob
import pandas
import pprint
from Bio            import Entrez
from Bio            import SeqIO
from Bio.Blast      import NCBIXML
#################################################
#               Entrez Settings                 #
#################################################
Entrez.email = "mfoster11@mgh.harvard.edu"

#################################################
#               Checkpoint Bools                #
#################################################
GET_REFERENCE_SEQS           = False
FIX_REFERENCE_SEQS           = False
WRITE_GENES                  = False
BUILD_BLAST_DB               = False
RUN_BLAST                    = False
PARSE_BLAST_RESULTS          = False
MERGE_RESULTS                = False
FINALIZE_AND_WRITE_METADATA  = False
#################################################
#              Function Definitions             #
#################################################

def parse_blast_results(blast_results : str) -> dict:
    file_size = os.path.getsize(blast_results)
    if file_size == 0:
        return None
    results = {}
    with open(blast_results, "r") as f:
        blast_records = NCBIXML.parse(f)
        for record in blast_records:
            results[record.query] = []
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    results[record.query].append((alignment.title, hsp.expect, hsp.score,))
    return results
#################################################
#         Define Inputs and File Paths          #
#################################################
assemblies = glob.glob("../../assemblies/*/annotation/*/*.gbff")
#print(assemblies)
metadata = "../metadata/longread_metadata_v2.csv"
input_metadata = load_input_metadata(metadata)
known_types = input_metadata[['Isolate', 'OspC_Type']]
###
# ******************************************
# * DELETED A LOT OF UNNECESSARY CODE HERE *
# ******************************************
###

#################################################
# ok now lets pull all of the ospCs from all assemblies and then blast them against our local blastdb
ospc_identifiers = [
    #"S2/P23", # This may not be right. I need to figure out why some assemblies lack OspC..
    "P12",
    "ospC",
    "OspC",
    "Ospc",
    "ospc",
    "outer surface protein C",
    "Outer surface protein C",
    "outer surface lipoprotein C",
    "Outer surface lipoprotein C",
    "outer surface protein c",
]

ospc_genes = {}
for file in assemblies:
    #print(file)
    sample_id = file.split("/")[-1].split(".")[0].replace("_200","" )
    ospc_genes[sample_id] = []
    with open(file, "r") as f:
        record = SeqIO.parse(f, "genbank")
        for rec in record:
            #print(rec.id)
            for feature in rec.features:
                if feature.type == "CDS":
                    product = feature.qualifiers['product']
                    if any(x in product[0] for x in ospc_identifiers):
                        #feature.qualifiers['product']
                        # get sequence
                        #print(feature.location)
                        # Now lets make a new SeqRecord object for the gene and append it to our list of OspC Genes!
                        seq = feature.extract(rec.seq)
                        gene = SeqIO.SeqRecord(seq, id=product[0], name=product[0], description=f"gene: {product[0]} from isolate {sample_id}", dbxrefs=None)
                        ospc_genes[sample_id].append(gene)

#################################################

if WRITE_GENES is True:
    # Okay now lets pull the ospC genes from all of the assemblies and write them to a file for blasting!
    for sample in ospc_genes:
        #print(sample)
        #print(ospc_genes[sample])
        with open(f"temp_in/{sample}.osps.fasta", "w") as f:
            SeqIO.write(ospc_genes[sample], f"temp_in/{sample}.osps.fasta", "fasta")
else:
    print("skipping writing ospC genes to file")

#################################################

if RUN_BLAST is True:
    # ok now we can blast each file against our local blastdb
    from Bio.Blast.Applications import NcbiblastnCommandline
    for sample in ospc_genes:
        print(f"blasting {sample}")
        cmd = NcbiblastnCommandline(query=f"temp_in/{sample}.osps.fasta", db="blastdb/ospc/ospc", out=f"temp_out/{sample}.osps_blast.xml", outfmt=5)
        print(cmd)
        stdout, stderr = cmd()
        print(stdout)
        print(stderr)
else:
    print("skipping blasting")
    
# ok now we can parse the blast results and then add the ospc type to the metadata
# we will use the function defined at the beginning to parse the blast results!

if PARSE_BLAST_RESULTS is True:
    checked_types_ospc = pandas.DataFrame(columns=["Isolate", "ospc_newtype"])
    blast_results_ospc = glob.glob("temp_out/*-osps_blast.xml")

    for result in blast_results_ospc:
        sample_id = result.split("/")[-1].split(".")[0].replace("-osps_blast", "")
        isolate = sample_id
        ospc_type = "unknown"
        parsed_results = parse_blast_results(result)
        if parsed_results is not None:
            for result in parsed_results:
                if len(parsed_results[result]) > 0:
                    best_hit = parsed_results[result][0][0]
                    coords = parsed_results[result][0][0].find("OspC_Type-")
                    ospc_type = parsed_results[result][0][0][coords+10:coords+14]
                    #print(isolate, ospc_type)
                    checked_types_ospc.loc[len(checked_types_ospc)] = [isolate, ospc_type]
                    break
        else:
            #print(isolate, ospc_type)
            checked_types_ospc.loc[len(checked_types_ospc)] = [isolate, "unknown"]
else:
    print("skipping parsing blast results")

#################################################

# This is to check how many have been duplicated to see what is messing up above...
print("Length of parsed results: ",len(checked_types_ospc))
duplicatesOC = checked_types_ospc[checked_types_ospc['Isolate'].duplicated()]
print("Duplicates OC : ",duplicatesOC)

# Visual Confirmation of Successful Parsing
#pprint.pprint(known_types)
#pprint.pprint(checked_types_ospc)
#pprint.pprint(checked_types_AT)

#################################################

    output_metadata = pandas.merge(input_metadata, merged_table, on='Isolate', how='left')
    output_metadata.to_csv(f"ospC_alleles_{version}.csv", index=False)
    print("metadata written!")
    
