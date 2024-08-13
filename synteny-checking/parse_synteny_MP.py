# multi_parse_homology.py

# Import some tools
import os
import sys
import csv
import glob
from Bio import SeqIO
from collections import defaultdict
# CLI stuff for MP
import argparse
from tqdm import tqdm
from pathlib import Path
import traceback
import tempfile
import pickle
from concurrent.futures import ProcessPoolExecutor, as_completed

def parse_genbank(path):
    records = defaultdict()
    for record in SeqIO.parse(path, 'genbank'):
        records[record.id] = record
    return records

def get_alignment_files(alignments_dir):
    alignments = glob.glob(f'{alignments_dir}/**/align_coords.tsv', recursive=True)
    return alignments

def check_alignment(path):
    # Let's make sure there's actually alignments within this file. Many such cases of no alignment. (Expected)
    # {{TO-DO: Add way to gather no-alns into a table}}
    with open(path, 'r') as infile:
        lines = infile.readlines()
    if len(lines) == 1:
        return False
    else:
        return True

def parse_alignment(path):
    # Read the file and parse the table.
    alignments = []
    with open(path, 'r') as infile:
        lines = infile.readlines()
    keys = lines[0].strip().split('\t')
    for line in lines[1::]:
        values = line.strip().split('\t')
        alignments.append({key:value for key,value in zip(keys, values)}) #zippity split
    return alignments

def is_within_range(feature, start, end):
    #print(feature)
    start_check = int(feature.location.start) >= start
    end_check = int(feature.location.end) <= end
    checked = (start_check and end_check) # gotta be within the range! {{TO-DO: Implement partial gene hit identification}}
    return checked

def get_features_from_range(record, start, end):
    # gonna pull the whole feature
    features = [feature for feature in record.features if is_within_range(feature, start, end) and feature.type == 'CDS']
    return features

def simplify_genes_for_contig(features):
    genes = []
    for feature in features:
        locus_tag = feature.qualifiers['locus_tag'][0].strip("'").strip('[').strip(']') # thanks python.
        product = ' '.join(feature.qualifiers['product'])
        gene = (locus_tag, product)
        genes.append(gene)
    return genes

def get_coverage(length, coords):
    covered_positions = set()
    for position in coords:
        start = position[0]
        end = position[1]
        if start > end:  # Note, this should only be required for calculation of coverage, feature extraction shouldn't need it.
            start, end = end, start
        covered_positions.update(range(start, end+1))
    percent_coverage = (len(covered_positions)/length)*100
    return percent_coverage

def get_homologies(alignments, asm1_dict, asm2_dict):
    # get_genes_from_alignments() => get_homologies()
    genes_dict = defaultdict() # set up dict for alignments
    for alignment in alignments: # iterate through each alignment for this particular comparison.
        # for whatever reason, mummer indicates the reference as the query and the assembly as the ref. Whatever. just be mindful.
        # Flipped for b31 alns, {{TO-DO: FIX THIS!}}
        # FLIPPED AGAIN FOR AVA >:(

        asm1_id        = alignment['QUERY_ID']
        asm1_name      = alignment['QUERY_NAME']
        asm1_start           = int(alignment['QUERY_START'])
        asm1_end             = int(alignment['QUERY_END'])
        asm1_aln_length      = int(alignment['QUERY_LENGTH']) # ALIGNED LENGTH)
        asm1_length    = len(asm1_dict[asm1_name].seq)
        #asm2_id             = alignment['REF_ID']
        asm2_name      = f'{alignment['REF_ID']}_{alignment['REF_NAME']}'#.replace("0000", "_").replace("_0", "_")
        asm2_start           = int(alignment['REF_START'])
        asm2_end             = int(alignment['REF_END'])
        asm2_aln_length      = int(alignment['REF_LENGTH']) # ALIGNED LENGTH)
        asm2_identity  = alignment['IDENTITY']
        asm2_length    = len(asm2_dict[asm2_name].seq)
        asm1_features = get_features_from_range(asm1_dict[asm1_name], asm1_start, asm1_end)
        asm2_features = get_features_from_range(asm2_dict[asm2_name], asm2_start, asm2_end)
        asm1_genes    = simplify_genes_for_contig(asm1_features)
        asm2_genes    = simplify_genes_for_contig(asm2_features)
        #asm_genes = "placeholder :)"
        # Get percent coverage for each alignment
        # UPDATE DON'T DO THAT HERE THE RANGES GET ALL WEIRD.
        #ref_coverage = get_coverage(ref_length, ref_start, ref_end)
        #contig_coverage = get_coverage(contig_length, contig_start, contig_end)

        if asm1_name not in genes_dict: # if the contig is not already in the dict, add it, also add contig_len to its own key.
            genes_dict[asm1_name] = defaultdict(dict)
            genes_dict[asm1_name]['contig_length'] = int(asm1_length) # Honestly for simplicity I should prob just make a separate dict for this.
            genes_dict[asm1_name]['assembly_id'] = asm1_id

        if asm2_name not in genes_dict[asm1_name]: # if the ref is not in the contig dict, add it and set val to an empty list.
            genes_dict[asm1_name][asm2_name] = []

        alignment_dict = {
                        'asm_aln': {
                                'start': int(asm1_start),
                                'end': int(asm1_end),
                                'aln_length': asm1_aln_length,
                                'percent_cov': int(),
                                'features': asm1_genes,
                            },
                        'ref_aln': {
                                'start': int(asm2_start),
                                'end': int(asm2_end),
                                'aln_length': int(asm2_aln_length),
                                'aln_identity': asm2_identity, # to the reference, this will probably get confusing downstream :)
                                'percent_cov': int(),
                                'features': asm2_genes,
                                'ref_length': int(asm2_length),
                            },
                    } # def the dictionary to append to the particular ref alignment for this single contig.

        genes_dict[asm1_name][asm2_name].append(alignment_dict)
        #print(genes_dict.keys())
    return genes_dict

def make_table_for_asm(genes_dict):
    lines = []
    for asm1_id in genes_dict.keys():
        asm1_genes = []
        asm2_genes = []
        for asm2_id in genes_dict[asm1_id].keys():
            assembly_id = genes_dict[asm1_id]['assembly_id']
            asm1_region_genes = []
            asm2_region_genes = []
            if asm2_id != 'contig_length' and asm2_id != 'assembly_id': # see above to-do re: separate dict just for len/assembly stats.
                asm1_len = genes_dict[asm1_id]['contig_length']
                asm2_coords = [(cov['asm_aln']['start'],cov['asm_aln']['end']) for cov in genes_dict[asm1_id][asm2_id]]
                asm1_aln_coords = [(cov['ref_aln']['start'],cov['ref_aln']['end']) for cov in genes_dict[asm1_id][asm2_id]]
                asm2_len = genes_dict[asm1_id][asm2_id][0]['ref_aln']['ref_length']
                for index, item in enumerate(genes_dict[asm1_id][asm2_id]):
                    asm1_region_genes = [gene[0] for gene in item['ref_aln']['features']]
                    asm2_region_genes = [gene[0] for gene in item['asm_aln']['features']]
                    asm1_genes.append(asm1_region_genes)
                    asm2_genes.append(asm2_region_genes)
                asm2_cov = get_coverage(asm2_len, asm2_coords)
                #print(contig_cov, aln_coords)
                asm1_cov = get_coverage(asm1_len, asm1_aln_coords)
                output_row = [assembly_id, asm1_id,asm1_len,asm2_id,asm2_len,f'{asm1_cov:.2f}',f'{asm2_cov:.2f}',len(asm1_genes),len(asm2_genes)]
                lines.append(output_row)
        return lines

def get_rows(alignments, genes_dict):
    """ Take the alns and format rows for dumping to a table separate from the main output"""
    rows = []
    for contig in genes_dict.keys():
        contig_list = []
        contig_len = genes_dict[contig]['contig_length']
        assembly_id = genes_dict[contig]['assembly_id']
        total_contig_cov = 0
        for aln in genes_dict[contig].keys():
            if aln != 'contig_length' and aln != 'assembly_id':
                asm_coords = [(cov['asm_aln']['start'],cov['asm_aln']['end']) for cov in genes_dict[contig][aln]]
                ref_coords = [(cov['ref_aln']['start'],cov['ref_aln']['end']) for cov in genes_dict[contig][aln]]
                ref_len = genes_dict[contig][aln][0]['ref_aln']['ref_length']
                contig_cov = get_coverage(contig_len, asm_coords)
                #collapsed_coords = collapse_coords(asm_coords)
                #asm_location = check_location_of_alignment(contig_len, asm_coords)
                #print(check_location_of_alignment(contig_len, collapse_coords(asm_coords)))
                ref_cov = get_coverage(ref_len, ref_coords)
                ref_aln = (f'({aln} : {contig_cov:.2f})')
                contig_list.append((ref_aln))
                total_contig_cov += contig_cov
        row = [assembly_id, contig ,contig_len, f'{total_contig_cov:.2f}']
        row.append(contig_list)
        rows.append(row)
    return rows

def write_rows(rows, output_file):
    with open(output_file, 'w') as outfile:
        writer = csv.writer(outfile, delimiter='\t')
        writer.writerow(['assembly_id', 'contig_id','contig_len','total_contig_coverage','list_of_alignments(ref:contig_cov:location)'])
        for row in rows:
            writer.writerow(row)

def parse_b31_alns(alignment, annotations_dir, ref_dict, output_dir):
    sample_id = alignment.split('/')[1]
    #print(f'parsing alignments for {sample_id}!')
    assembly = f'{annotations_dir}/{sample_id}/{sample_id}.gbff'
    print(f'gbff for {sample_id} exists: {os.path.exists(assembly)}')
    #print(assembly)
    asm_dict = parse_genbank(assembly)
    alignments = parse_alignment(alignment)
    homologies = get_homologies(alignments, asm_dict, ref_dict)
    output_file_detailed = f'{output_dir}/detailed_coverage/{sample_id}_B31_Synteny.tsv'
    make_table_for_asm(homologies, output_file_detailed)
    rows = get_rows(alignments, homologies)
    detailed_out = f'{output_dir}/detailed_coverage/'
    if not os.path.exists(detailed_out):
        os.makedirs(detailed_out)
    simple_out = f'{output_dir}/simple_coverage/'
    if not os.path.exists(simple_out):
        os.makedirs(simple_out)
    output_file_simple = f'{output_dir}/simple_coverage/{sample_id}_coverage.tsv'
    write_rows(rows, output_file_simple)

def run_all_b31_alns(alignments_dir, annotations_dir, reference_genbank, output_dir):
    print(f'annotations dir: {annotations_dir}')
    print('parsing ref genbank into dictionary')
    ref_dict = parse_genbank(reference_genbank)
    alignment_files = glob.glob(f'{alignments_dir}/*/align_coords.tsv', recursive=True)
    print(f'found {len(alignment_files)} alignment files')
    for file in alignment_files:
        print(f'parsing {file}')
        parse_b31_alns(file, annotations_dir, ref_dict, output_dir)
        print(f'Finished! moving on')

def parse_ids_from_filename(alignment):
    alignment_dir_name = os.path.dirname(alignment).split('/')[-1]
    asm1_id = os.path.basename(alignment_dir_name).split('_vs_')[0]
    asm2_id = os.path.basename(alignment_dir_name).split('_vs_')[1]
    return asm1_id, asm2_id

def parse_single_pair_aln(alignment_file, annotations_dir):
    # {{TO-DO: Not at this point but at some point I need to pull in the plasmid ID to name mapping dict.}}
    ## okay let's get our ids.
    asm1_id, asm2_id = parse_ids_from_filename(alignment_file)

    ## First let's check to see that there are actually alignments for this pair.
    completion_msg = f'Parsed homology between {asm1_id} and {asm2_id}!'
    if not check_alignment(alignment_file):
        completion_msg = f'NO HOMOLOGY BETWEEN {asm1_id} and {asm2_id}!'
        return (False, completion_msg, [asm1_id, asm2_id]) # return the two ids if no homology!

    ## now let's parse the genbanks
    asm1_gb = f'{annotations_dir}/{asm1_id}.gbff'
    asm2_gb = f'{annotations_dir}/{asm2_id}.gbff'
    asm1_dict = parse_genbank(asm1_gb)
    asm2_dict = parse_genbank(asm2_gb)

    ### ok, it's not empty, let's parse this out.
    alignments = parse_alignment(alignment_file)

    ## Anyway let's get our genes for these alignments and actually *parse* the alignments.
    homologies = get_homologies(alignments, asm1_dict, asm2_dict)
    detailed_rows = make_table_for_asm(homologies)
    ## Now let's make our rows.
    simple_rows = get_rows(alignments, homologies)
    completion_msg = f'Parsed homology between {asm1_id} and {asm2_id}!'

    return (True, completion_msg, simple_rows, detailed_rows)

def catch_error(function, *args, **kwargs):
    try:
        results = function(*args, **kwargs)
        return results
    except Exception as e:
        error_message = f'{str(e)}\n{traceback.format_exc()}'
        return error_message


def parallel_parse(cpus, alignment_files, annotations_dir):
    # def parse_all_v_all():
    # this may require parallelization?
    # yeah let's just go ahead and do that.

    simple_rows = []
    detailed_rows = []
    no_homology = []

    no_homology_header = ['ASM1', 'ASM2']
    detailed_header_row = ['contig_id', 'contig_len', 'ref', 'ref_len', 'contig_cov', 'reference_cov', 'genes_on_ref', 'genes_on_contig']
    simple_header = ['contig_id','contig_len','total_contig_coverage','list_of_alignments(ref:contig_cov:location)']

    no_homology_header.append([no_homology_header])
    simple_rows.append(simple_header)
    detailed_rows.append(detailed_header_row)

    with ProcessPoolExecutor(max_workers=cpus) as executor: # need to specify this elsewhere.
        futures = []
        # okay so we need to get all of our individual alignments, build a list of args, then feed em into the workers.
        # first let's parse the alignments and figure out how to divide this.
        # Wait, it's literally just iteration through a list.
        # anyway let's set up the big list o' rows.
        for alignment in alignment_files:
            futures.append(executor.submit(catch_error, parse_single_pair_aln, alignment, annotations_dir)) # feed our single command the aln and the dir for the genomes
        # ok now let's gather each alignment and cat it onto the list of rows!
        with tqdm(total=len(futures)) as pbar:
            for future in as_completed(futures):
                try:
                    result = future.result()
                    if result[0] is True:
                        simple_rows.extend(result[2])
                        detailed_rows.extend(result[3])
                    elif result[0] is False:
                        no_homology.extend([result[2]])
                    else:
                        print(f"ERROR WITH PARSING: {result}")
                except Exception as e:
                    tqdm.write(f"Error: {e}")# double it so we can keep a record on the screen.
                pbar.update(1)
                #pbar.write(result[1])

    return simple_rows, detailed_rows, no_homology

def custom_write(text):
    # Use the tqdm.write method to ensure that the progress bar does not get disrupted
    tqdm.write(text)
    # Move the cursor up one line and clear the line
    sys.stdout.write('\033[F\033[K')

def search_in_subdir(directory):
    return list(directory.rglob("align_coords.tsv"))

def parallel_find(directory, cpus):
    # get list of results directories.
    dir_path = Path(directory)
    files = []
    futures = []
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        subdirs = [subdir for subdir in dir_path.iterdir() if subdir.is_dir()]
        futures = [executor.submit(search_in_subdir, subdir) for subdir in subdirs]
        for future in tqdm(as_completed(futures), total=len(futures), desc="finding results files", unit=" files "):
            result = future.result()
            files.extend(result)
    return files


def main():
    # Create the parser
    parser = argparse.ArgumentParser(description="A script to run mauve on a directory of assemblies against the B31 reference genome.")
    # Add the arguments
    parser.add_argument('alignments_dir',  help='The directory containing the alignments to parse', type=str)
    parser.add_argument('annotations_dir', help='The directory containing the annotations (genbanks)', type=str)
    parser.add_argument('output_dir',      help='The directory for outputs', type=str)
    parser.add_argument('--cpus',          help='How many cores we rippin', type=int)
    parser.add_argument('--mode',          help="specify \'ava\' for all_vs_all and 'av1' for all_vs_one mode, 'av1' requires you to give the reference genbank with --reference", type=str)
    parser.add_argument('--ref',           help='the reference genbank for all_vs_one mode', required=False, type=str)
    parser.add_argument('--resume',        help='use this to resume without searching for all of the results again.', required=False, action='store_true',)

    # Parse the arguments
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    mode = args.mode
    files_list_name = 'results_file_paths.pkl'
    pickle_simple_rows = 'parsed_simple_rows.pkl'
    pickle_detailed_rows = 'parsed_detailed_rows.pkl'
    pickle_no_homology = 'parsed_no_homology.pkl'
    files_list_path = os.path.join(args.output_dir, files_list_name)
    pickle_simple_rows = os.path.join(args.output_dir, pickle_simple_rows)
    pickle_detailed_rows = os.path.join(args.output_dir, pickle_detailed_rows)
    pickle_no_homology =  os.path.join(args.output_dir, pickle_no_homology)

    if args.resume:
        print(f"Reading list of results files from {files_list_path}")
        alignment_files = pickle.load(open(files_list_path,'rb'))
    else:
        alignment_files = parallel_find(args.alignments_dir, args.cpus) # this is so massive it needed parallelization just to build the list of files :(
        print(f"writing list of results files to  {files_list_path}")
        pickle.dump(alignment_files, open(files_list_path, 'wb'))

    #alignment_files = get_alignment_files(args.alignments_dir)
    if mode == 'ava':
        if args.resume:
            print(f"Unpacking parsed pickles for processing:\n{pickle_simple_rows}\n{pickle_detailed_rows}\n{pickle_no_homology}")
            simple_rows   = pickle.load(open(pickle_simple_rows, 'rb'))
            detailed_rows = pickle.load(open(pickle_detailed_rows, 'rb'))
            no_homologies = pickle.load(open(pickle_no_homology, 'rb'))
            print("Pickles Popped!")
        else:
            simple_rows, detailed_rows, no_homologies = parallel_parse(args.cpus, alignment_files, args.annotations_dir)
            print(f"Promptly packing pickles of parsed products:\n{pickle_simple_rows}\n{pickle_detailed_rows}\n{pickle_no_homology}")
            pickle.dump(simple_rows, open(pickle_simple_rows, 'wb'))
            pickle.dump(detailed_rows, open(pickle_detailed_rows, 'wb'))
            pickle.dump(no_homologies, open(pickle_no_homology, 'wb'))
            print("Pickles Packed!")
            ## Now let's make our rows.
        #simple_rows = get_rows(alignments, homologies)
        #completion_msg = f'Parsed homology between {asm1_id} and {asm2_id}!'
        output_simple = f'{args.output_dir}/ava_homo_simple.tsv'
        output_detailed = f'{args.output_dir}/ava_homo_detailed.tsv'
        output_no_homo = f'{args.output_dir}/ava_no_homo.tsv'
        #print(detailed_rows[1:5])
        write_rows(simple_rows,   output_simple)
        write_rows(detailed_rows, output_detailed)
        write_rows(no_homologies, output_no_homo)
        return print('Finished!')

    if mode == 'av1':
        print("running in av1 mode")
        run_all_b31_alns(args.alignments_dir, args.annotations_dir, args.ref, args.output_dir)
        return 0

if __name__ == "__main__":
    main()