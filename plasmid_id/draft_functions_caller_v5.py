## Parsing stuff

#def get_parsing_table(pickle_path):
#    return pickle.load(open(f'{parsing_tables_dir}/blast_parsing_dict.pkl', 'rb'))

#def get_hit_name(alignment_id):
#    """Extract hit id from hit title"""
#    if "|" in alignment_id:
#        if 'pdb' in alignment_id:
#            print(alignment_id)
#            plasmid_name = alignment_id.split("|")[-1]
#            print(plasmid_name)
#            return plasmid_name
#        else:
#            plasmid_name = alignment_id.split("|")[1]
#            return plasmid_name
#    else:
#        plasmid_name = alignment_id
#        return plasmid_name
#    return plasmid_name
    
#def parse_results_file(xml_file, parsing_dict):
#    """ JUST FOR PF32 RESULTS FOR NOW!!!"""
#    parsed_hits = defaultdict(dict)
#    with open(xml_file, 'r') as file:
#        records = NCBIXML.parse(file)
#    # leterrip
#    for record in records:
#        query = record.query
#        if query not in parsed_hits:
#            parsed_hits[query] = []
#        hits = parsed_hits[query]
#        query_length = record.query_length
#        for alignment in records.alignments:
#            alignment_id = alignment.hit_id
#            # name parsing shenanigans, probably broken!
#            ref_name = get_hit_id(hit_id).split('_')[-1] if len(plasmid_name.split('_')) > 1 else plasmid_name
#            strain = alignment_id.split('_')[0]
#            # get ref length from our parsing dictionary!
#            for ncbi_id, data in parsing_table.items():
#                if data['strain'] == strain and data['name'] == ref_name:
#                    #print(ncbi_id)
#                    ref_length = parsing_table[ncbi_id]['length']
#                else:
#                    ref_length = 'NaN'
#            # Okay now lets get our coverage details :)
#            overall_percent_identity, coverage_percentage, covered_length, covered_intervals, query_intervals, subject_hit_coords = calculate_percent_identity_and_coverage(alignment)
#            alignment_dict = {
#                "alignment_id": alignment_id,
#                "plasmid_name": ref_name,
#                "query_length": query_length,
#                "ref_length": ref_length,
#                "ref_total_length": ref_total_length,
#                "overall_percent_identity": overall_percent_identity,
#                "coverage_percentage": coverage_percentage,
#                'covered_positions': covered_length,
#                'covered_intervals': covered_intervals,
#                'query_intervals': query_intervals,
#                'subject_hit_coords': subject_hit_coords,
#            }
#            hits.append(alignment_dict)
#    return parsed_hits # this is a dict with structure: { query1: [hit, hit, ...], query2: [...], ...}
#
#def write_pf32_hits(assembly_id, parsed_hits, outfile):
#    with open(outfile, 'w', newline='') as tsvfile:
#        fieldnames = ['assembly', 'contig', 'contig_len', 'pf32', 'pf32_ident', 'pf32_coverage']
#        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
#        writer.writeheader()
#
#        for contig, hits in parsed_hits.items():
#            pf32_best_identity = get_best_match(hits), 'overall_percent_identity')
#            pf32_best_coverage = get_best_match(hits), 'coverage_percentage')
#
#            best_pf = None
#            contig_len = ''
#            
def main():
    #### THIS IS A WEE BIT TOO COMPLICATED TO IMPLEMENT BEFORE EVEN GETTING THE THING UP AND RUNNING :(
    single_input = False
    parse = False
    blastem = False
    # Create the parser
    parser = argparse.ArgumentParser(prog = "PlasmidCaller_v5",
                                    description = "{{TO-DO: POPULATE DESCRIPTION!}}")
    # Add the mandatory arguments for all modes
    parser.add_argument('--cpus',        required=True, type=int,  help='How many cores we rippin')
    parser.add_argument('--input',       required=True, type=str,  help='A directory/single input file')
    parser.add_argument('--input_type',  required=True, type=str,  help='what is the type of the input file(s)? genbank or fasta')
    parser.add_argument('--annotations', required=True, action= ,type=str,  help='The directory/single file containing the annotations (genbanks)')
    parser.add_argument('--output',      required=False, action= ,type=str,  help='The directory for outputs') # ok this may not be initially required.
    # init subparsers for each mode of running this script
    subparsers = parser.add_subparsers(dest='command', help='sub-command help')

    # subparser for making a database
    make_db_parser = subparsers.add_parser(name='--make_db', help='Make a blast database using the given input')
    make_db_parser.add_argument('--db_type', type=str, help='type of blast db to create \'nucl\' or \'prot\'')
    #make_db_parser.add_argument('--db_dir',  type=str, help='directory to output the database after creation.') # this is covered by --output
    make_db_parser.add_argument('--db_name', type=str, help='desired name of the database.')

    # subparser for running blast.
    run_blast_parser = parser.add_subparsers(name='--blast', help='Analyze inputs against given databases using blast!')
    run_blast_parser.add_argument(name='--db', action='append', type=str, help='database to use in analysis, for multiple databases: --db path/to/db1 --db path/to/db2')

    # subparser for parsing blast results
    parse_results = parser.add_subparsers(name='--parse', help='Parse blast results and output to table!')
    parse_results.add_argument('--parsing_tables', help='directory for useful parsing tables') # this may not be useful? needs work! 



    # Parse the arguments
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)

    # Determine if single input or multiple inputs.
    if Path(args.input).suffix:
        # gotta be one of these allowed inputs.
        if Path(args.input).suffix not in acceptable_inputs:
            parser.error(f"Please specify one of the following input files: {acceptable_inputs}")
            return 0
        else:
            single_input = True # change flag to True

    if args.command == 'make_db':
        input_seqs = args.input
        db_type = args.db_type
        out_dir = args.output
        db_name = args.db_name
        db_out = Path.join(out_dir, db_name)
        command = make_blast_db(db_type, input_seqs, db_out)
        run_command(command)
        return 0
    elif args.command == 'blast':
        if single_input:
            seq = args.input
        else:
            seqs = get_input_files()    
     