#!/usr/bin/env python3
#################################################
from gzip import WRITE
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
GET_REFERENCE_SEQS = False
FIX_REFERENCE_SEQS = False
WRITE_GENES = False
BUILD_BLAST_DB = False
RUN_BLAST = False
PARSE_BLAST_RESULTS = False
MERGE_RESULTS = False
FINALIZE_AND_WRITE_METADATA = False
#################################################
#              Function Definitions             #
#################################################
def load_input_metadata(metadata : str) -> pandas.DataFrame:
    input_metadata = pandas.read_csv(metadata)
    return input_metadata

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
#print(known_types)
#################################################
# # Okay, Ira has sent me a large file of OspC types and AA seqs
# I will be taking this, making a blast DB and then I need to isolate the OspC seqs from each assembly, then blast against the OspC database.
# Then parse the blast results and add the type to the metadata table.
#################################################
# OspC List from Ira and Hanincova et al.
ospc_genotype_to_ref = {
     'A' : 'X69596',   'Ba' : 'EF537413',
    'Bb' : 'NC_011724', 'C' : 'DQ437462',
    'Da' : 'AF029863', 'Db' : 'GQ478283',
     'E' : 'AY275221', 'Fa' : 'AY275225',
    'Fb' : 'EF537433', 'Fc' : 'GQ478285',
     'G' : 'AY275223', 'Ha' : 'EU377781',
    'Hb' : 'GQ478286', 'Ia' : 'AY275219',
    'Ib' : 'EU377752',  'J' : 'CP001535',
     'K' : 'AY275214',  'L' : 'EU375832',
     'M' : 'CP001550',  'N' : 'EU377775',
     'O' : 'FJ997281',  'T' : 'AY275222',
    'Ua' : 'EU377769', 'Ub' : 'GQ478287',
    'A3' : 'EF592541', 'B3' : 'EF592542',
    'C3' : 'EF592543', 'D3' : 'EF592544',
    'E3' : 'EF592545', 'F3' : 'EF592547',
    'H3' : 'FJ932733', 'I3' : 'FJ932734',}

# OspC AT list from file provided by Ira
ospc_AT_to_ref = {
   'AT1' : 'EU482041', 'AT2' : 'EU482042',
   'AT3' : 'EU482043', 'AT4' : 'EU482044',
   'AT5' : 'EU482045', 'AT6' : 'EU482046',
   'AT7' : 'EU482047', 'AT8' : 'EU482048',
   'AT9' : 'EU482049', 'AT10' : 'EU482050',
  'AT11' : 'EU482051', 'AT12' : 'EU482052',
  'AT13' : 'EU482053', 'AT14' : 'EU482054',
  'AT15' : 'EU482055', 'AT16' : 'EU482056', 
  'B.bissettii_25015' : 'U04282'}

#################################################

if GET_REFERENCE_SEQS is True:
    # ok now to pull all of those seqs via entrez direct and then make a multifasta file upon which to build a blastdb from
    # we will use the following command to pull the sequences from the NCBI nucleotide database using our defined dictionary
    # and we will also call the command using the biopython wrapper for entrez direct
    for ospc_type, ref in ospc_genotype_to_ref.items():
        print(f"fetching record for {ospc_type} : {ref}")
        handle = Entrez.efetch(db="nucleotide", id=ref, rettype="fasta", retmode="text")
        record = handle.read()
        with open(f"ospc_seqs/{ospc_type}.fasta", "w") as f:
            print("writing record to file!")
            f.write(record[0::])
        handle.close()
    print("done!")
    for ospc_AT, ref in ospc_AT_to_ref.items():
        print(f"fetching record for {ospc_AT} : {ref}")
        handle = Entrez.efetch(db="nucleotide", id=ref, rettype="fasta", retmode="text")
        record = handle.read()
        with open(f"AT_seqs/{ospc_AT}.fasta", "w") as f:
            print("writing record to file!")
            f.write(record[0::])
        handle.close()
    print("done!")
else:
    print("skipping reference sequence retrieval")

#################################################

if FIX_REFERENCE_SEQS is True:
    # ok so some of those entries include the entire plasmid sequence not just the ospc gene, so we will need to extract the ospc gene from the plasmid sequence
    # the ones that are plasmid sequences are as follows:
    # Bb, J, M,
    # we will need to extract the ospc gene from these sequences
    for ospc_type in ['Bb', 'J', 'M']:
        print(ospc_genotype_to_ref[ospc_type])

    coords = {
        'Bb' : (16904,17540), # DON'T FORGET TO CHANGE TO 0 BASED COORDINATES!
        'J' : (16909,17545), # DON'T FORGET TO CHANGE TO 0 BASED COORDINATES!
        'M' : (16916,17555), # DON'T FORGET TO CHANGE TO 0 BASED COORDINATES!
    }
    #
    ## ok we already have the full plasmid seqs for each of these three types, we also now have the coordinates
    for ospc_type in ['Bb', 'J', 'M']:
        print(f"extracting ospc gene from {ospc_type} sequence")
        with open(f"ospc_seqs/{ospc_type}.fasta", "r") as f:
            record = SeqIO.read(f, "fasta")
            print(record)
            ospc_gene = record[coords[ospc_type][0]:coords[ospc_type][1]]
            with open(f"ospc_seqs/{ospc_type}_ospc.fasta", "w") as f:
                SeqIO.write(ospc_gene, f, "fasta")
else:
    print("skipping ospc gene extraction and record repair")

#################################################

if WRITE_REF_SEQS is True:
    # Now lets write out all of our OspC sequences to a single file for blastdb creation
    ospc_seqs = glob.glob("ospc_seqs/*.fasta")
    #print(ospc_seqs)
    ospc_seq_file = []
    for seq in ospc_seqs:
        record = SeqIO.read(seq, "fasta")
        record.id = "OspC_Type-" + seq.split("/")[-1].split(".")[0]
        record.description = ""
        #print(record.id)
        ospc_seq_file.append(record)
    with open("all_ospc.fasta", "w") as f:
        SeqIO.write(ospc_seq_file, f, "fasta")
    print("OspC file written!")
else:
    print("skipping writing out reference OspC sequences")

if WRITE_REF_SEQS is True:
    # Okay now lets write the AT sequences to a single file for blastdb creation
    AT_seqs = glob.glob("AT_seqs/*.fasta")
    #print(AT_seqs)
    AT_seq_file = []
    for seq in AT_seqs:
        record = SeqIO.read(seq, "fasta")
        record.id = "OspC_Type-" + seq.split("/")[-1].split(".")[0]
        record.description = ""
        #print(record.id)
        AT_seq_file.append(record)
    with open("all_ospc_ATs.fasta", "w") as f:
        SeqIO.write(AT_seq_file, f, "fasta")
    print("AT file written!")
else:
    print("skipping writing out reference AT sequences")

#################################################

if MAKE_BLASTDB is True:
    # ok now lets make a local blast DB for our OspC types
    from Bio.Blast.Applications import NcbimakeblastdbCommandline
    cmd = NcbimakeblastdbCommandline(input_file="all_ospc.fasta", dbtype="nucl", out="blastdb/ospc/ospc")
    print(cmd)
    stdout, stderr = cmd()
    print(stdout)
    print(stderr)
    # now lets make a local blast DB for our AT types
    from Bio.Blast.Applications import NcbimakeblastdbCommandline
    cmd = NcbimakeblastdbCommandline(input_file="all_ospc_ATs.fasta", dbtype="nucl", out="blastdb/ospcAT/ospcAT")
    print(cmd)
    stdout, stderr = cmd()
    print(stdout)
    print(stderr)
else:
    print("skipping making blastdb for AT")

#################################################
# ok now lets pull all of the ospCs from all assemblies and then blast them against our local blastdb
ospc_identifiers = [
    #"S2/P23", # This may not be right. I need to figure out why some assemblies lack OspC.....
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

if RUN_BLAST is True:
    # ok now we can blast each file against our local blastdb
    from Bio.Blast.Applications import NcbiblastnCommandline
    for sample in ospc_genes:
        print(f"blasting {sample}")
        cmd = NcbiblastnCommandline(query=f"temp_in/{sample}.osps.fasta", db="blastdb/ospcAT/ospcAT", out=f"temp_out/{sample}.ospsAT_blast.xml", outfmt=5)
        print(cmd)
        stdout, stderr = cmd()
        print(stdout)
        print(stderr)
else:
    print("skipping blasting")

#################################################

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
    print("skipping parsing blast OC results")

if PARSE_BLAST_RESULTS is True:
    checked_types_AT = pandas.DataFrame(columns=["Isolate", "AT_type"])
    blast_results_AT = glob.glob("temp_out/*-ospsAT_blast.xml")

    for result in blast_results_AT:
        #print(result)
        sample_id = result.split("/")[-1].replace("-ospsAT_blast.xml", "")
        isolate = sample_id
        parsed_results = parse_blast_results(result)
        if parsed_results is not None:
            for result in parsed_results:
                if len(parsed_results[result]) > 0:
                    best_hit = parsed_results[result][0][0]
                    coords = parsed_results[result][0][0].find("OspC_Type-")
                    ospc_type = parsed_results[result][0][0][coords+10:coords+14]
                    #print(isolate, ospc_type)
                    checked_types_AT.loc[len(checked_types_AT)] = [isolate, ospc_type]
                    break
        else:
            #print(isolate, ospc_type)
            checked_types_AT.loc[len(checked_types_AT)] = [isolate, "unknown"]
else:
    print("skipping parsing blast AT results")

#################################################

# This is to check how many have been duplicated to see what is messing up above...
print("Length of parsed results: ",len(checked_types_ospc))
print("Length of parsed results: ",len(checked_types_AT))
duplicatesOC = checked_types_ospc[checked_types_ospc['Isolate'].duplicated()]
print("Duplicates OC : ",duplicatesOC)
duplicatesAT = checked_types_AT[checked_types_AT['Isolate'].duplicated()]
print("Duplicates AT : ",duplicatesAT)

# Visual Confirmation of Successful Parsing
#pprint.pprint(known_types)
#pprint.pprint(checked_types_ospc)
#pprint.pprint(checked_types_AT)

#################################################

if MERGE_RESULTS is True:
    merged_table = []
    # Merge the two tables on the 'Isolate' column
    merged_table = pandas.merge(known_types, checked_types_ospc, on='Isolate', how='left')
    #print(merged_table)
    merged_table = pandas.merge(merged_table, checked_types_AT, on='Isolate', how='left')
    #print(merged_table)
    for i in range(0, len(merged_table)):
        if pandas.isnull(merged_table.at[i, 'OspC_Type']):
            merged_table.at[i, 'OspC_Type'] = merged_table.at[i, 'ospc_newtype']
    merged_table.drop(columns=['ospc_newtype'], inplace=True)
    #print(merged_table)
    print("Known types and new types tables have been merged!")
else:
    print("skipping merging")

#################################################

if FINALIZE_AND_WRITE_METADATA is True:
    #input_metadata.drop(columns=['OspC_Type'], inplace=True)
    output_metadata = pandas.merge(input_metadata, merged_table, on='Isolate', how='left')
    #print(output_metadata)
    output_metadata.to_csv("longread_metadata_v3.csv", index=False)
    print("metadata written!")
else:
    print("skipping writing final metadata")
