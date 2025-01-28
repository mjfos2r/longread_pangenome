import os
import glob
import pandas
import tempfile
from pathlib import Path
from Bio import SeqIO
from collections import defaultdict
from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.align import AlignCoord
from pygenomeviz.utils import ColorCycler

# Lets start with some functions to process the contig, this will come from iterating through the align_coords.tsv dataframe
# and grouping by contig_id.

def process_assembly(homologies_file, syntenies_file, reference_gb, assembly_gb, output_dir, **kwargs):
    # process the alignments and spit out a plot for each contig :)
    gbk_assembly = list(SeqIO.parse(open(assembly_gb, 'r'), 'genbank'))
    gbk_reference = list(SeqIO.parse(open(reference_gb, 'r'), 'genbank'))
    assembly_id = Path(assembly_gb).stem
    reference_id = '.'.join(os.path.basename(reference_gb).split('.')[:-1])

    output_path = os.path.join(output_dir, assembly_id)
    png_path = os.path.join(output_path, 'png')
    html_path = os.path.join(output_path, 'html')

    check_n_make(output_path)
    check_n_make(png_path)
    check_n_make(html_path)

    homologies = pandas.read_csv(homologies_file, sep='\t')
    syntenies = pandas.read_csv(syntenies_file, sep = '\t')

    grouped_homologies = homologies.groupby('REF_NAME', ) # THIS WILL REQUIRE CHANGING IN THE FUTURE {{TO-DO: FIX FOR NON B31 CASE}}
    grouped_syntenies  = syntenies.groupby('REF_NAME', ) # THIS WILL REQUIRE CHANGING IN THE FUTURE {{TO-DO: FIX FOR NON B31 CASE}}

    for contig_id, contig_df in grouped_homologies:
        num_homologies = len(contig_df)
        if check_group_exists(grouped_syntenies, contig_id):
            num_syntenies = len(grouped_syntenies.get_group(contig_id))
        else:
            num_syntenies = 0
        print(f'Processing {contig_id}: homologies:{num_homologies}, syntenies:{num_syntenies}')
        if num_syntenies > 0:
            process_and_plot(contig_id, contig_df, grouped_syntenies, gbk_assembly, gbk_reference, assembly_id, reference_id, output_path)
        else:
            print(f"\tError plotting {contig_id}! No syntenies! Attempting to only plot homology!\n")
            process_and_plot(contig_id, contig_df, grouped_syntenies, gbk_assembly, gbk_reference, assembly_id, reference_id, output_path, homology_only=True)
    print(f"Finished plotting contigs for {assembly_id}! Have a wonderful day :)")

def process_and_plot(contig_id, contig_df, grouped_syntenies, gbk_assembly, gbk_reference, assembly_id, reference_id, output_path, **kwargs):
    # Process contig alignments for homology and synteny
    homology_only = kwargs.get('homology_only', False)
    # ok now let's parse this out,
    ref_homology_ids, homology_coords = process_contig_alignments(contig_df)
    if homology_only:
        reference_ids = ref_homology_ids
    else:
        ref_synteny_ids, synteny_coords = process_contig_alignments(grouped_syntenies.get_group(contig_id))
        reference_ids = ref_homology_ids.union(ref_synteny_ids)

    # Parse GenBank files for assembly and reference
    gbks = []
    gbk_asm = parse_split_genbank(gbk_assembly, contig_id, assembly_id)
    gbks.append(gbk_asm)
    gbk_ref = parse_split_genbank(gbk_reference, reference_ids, reference_id)
    gbks.append(gbk_ref)

    # Initialize GenomeViz
    gv = GenomeViz(track_align_type="center")

    # Add tracks for each GenBank file
    track_number = 0
    for gbk in gbks:
        if track_number == 0:
            sublabel_pos = "top-center"
        else:
            sublabel_pos = "bottom-center"
        color = ColorCycler()
        #print(gbk.name)
        track = gv.add_feature_track(gbk.name, gbk.get_seqid2size(), space=0.05, label_kws=dict(color=color), align_label=False)
        for seqid, features in gbk.get_seqid2features(feature_type="CDS").items():
            segment = track.get_segment(seqid)

            segment.add_sublabel(f"{segment.name}: {segment.start:,}-{segment.end:,}bp", size = 8, pos = sublabel_pos, ymargin=0.3)
            for feature in features:
                segment.add_features(feature, fc="blue", lw=0.05, ignore_outside_range = True, plotstyle="rbox")

        track_number += 1
    # Add homology and synteny links
    if len(homology_coords) > 0:
        min_ident = int(min([ac.identity for ac in homology_coords if ac.identity]))
        color, inverted_color = "grey", "red"
        for ac in homology_coords:
            gv.add_link(ac.query_link, ac.ref_link, color=color, inverted_color=inverted_color, filter_length=250, curve=True, v=ac.identity, vmin=min_ident)
        gv.set_colorbar([color, inverted_color], vmin=min_ident, bar_label="Identity")
    if not homology_only:
        color, inverted_color = "lightblue", "green"
        for sc in synteny_coords:
            gv.add_link(sc.query_link, sc.ref_link, color=color, inverted_color=inverted_color, curve=True)

    # Save the plot as PNG and HTML
    output_png = os.path.join(output_path, 'png', f'{contig_id}_synteny_plot.png')
    output_html = os.path.join(output_path, 'html', f'{contig_id}_synteny_plot.html')
    gv.savefig(output_png)
    gv.savefig_html(output_html)
    print(f'Finished plotting {contig_id}\n')

def process_contig_alignments(contig_df):
    #keys = ['QUERY_ID', 'QUERY_NAME', 'QUERY_START', 'QUERY_END', 'QUERY_LENGTH',
    #   'REF_ID', 'REF_NAME', 'REF_START', 'REF_END', 'REF_LENGTH', 'IDENTITY',
    #   'EVALUE']
    reference_ids = set()
    acs = []
    for index, row in contig_df.iterrows():
        row_dict = {k:v for k,v in row.items()}
        reference_ids.add(row_dict['QUERY_NAME']) # change this to query_id for all_v_all ? doesn't matter for B31
        acs.append(
            AlignCoord(
                row_dict['QUERY_ID'],
                row_dict['QUERY_NAME'],
                row_dict['QUERY_START'],
                row_dict['QUERY_END'],
                row_dict['REF_ID'],
                row_dict['REF_NAME'],
                row_dict['REF_START'],
                row_dict['REF_END'],
                row_dict['IDENTITY'],
            )
        )
    return reference_ids, acs

def parse_genbanks_for_contig(reference_ids):
    files = []
    for genbank_id in reference_ids:
        genbank_file = get_genbank_file(genbank_id)
        files.append(genbank_file)
    combined_temp_file = concatenate_genbank_files(files)
    gbk_ref = Genbank(combined_temp_file)
    return gbk_ref

def parse_split_genbank(records, record_ids, temp_file_name):
    temp_dir = tempfile.gettempdir()
    temp_file_path = os.path.join(temp_dir, temp_file_name)
    filtered_records = [record for record in records if record.id in record_ids]
    with open(temp_file_path, 'w') as outfile:
        SeqIO.write(filtered_records, outfile, 'genbank')
    gbk_ref = Genbank(temp_file_path, name=temp_file_name)
    return gbk_ref

def check_group_exists(grouped_df, group_name):
    return group_name in grouped_df.groups

def check_n_make(path):
    if not os.path.exists(path):
        os.makedirs(path)

def concatenate_genbank_files(file_list):
    # Create a temporary file
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.gb')

    with open(temp_file.name, 'w') as outfile:
        for file_name in file_list:
            with open(file_name, 'r') as infile:
                records = SeqIO.parse(infile, 'genbank')
                SeqIO.write(records, outfile, 'genbank')
    # Return the name of the temporary file
    return temp_file.name

def get_genbank_file(genbank_id, where_to_look):
    file_path = os.path.join(where_to_look, f'{genbank_id}.gbff')
    if os.path.exists(file_path):
        return genbank_file
    else:
        raise FileNotFoundError(f"The file at path '{file_path}' does not exist.")

def make_track_from_gb(gv, track_id, lens, features):
    track = gv.add_feature_track(track_id, lens)
    track.add_sublabel()
    track.add_features(features)
    return track

def main():

    #
    #
    #
    #
    # DO TOMORROW :)