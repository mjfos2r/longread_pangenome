import os
import sys
import csv
import glob
import time
from Bio import SeqIO
from collections import defaultdict
import argparse
from pathlib import Path
import traceback
import pickle
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeRemainingColumn, TimeElapsedColumn
from rich.console import Console
from rich.panel import Panel
from rich.live import Live
from rich.table import Table
from rich.layout import Layout
from rich import box

class ConsolePanel(Console):
    def __init__(self,*args,**kwargs):
        console_file = open(os.devnull,'w')
        super().__init__(record=True,file=console_file,*args,**kwargs)

    def __rich_console__(self,console,options,stderr=False):
        texts = self.export_text(clear=False,).split('\n') # styles=True <- This fucks up the borders.
        for line in texts[-options.height:]:
            yield line
        
console = ConsolePanel()

def parse_genbank(path):
    try:
        records = list(SeqIO.parse(path, 'genbank'))
        if not records:
            console.print(f"[bold yellow]Warning:[/bold yellow] Empty GenBank file: {path}")
            return None
        return records[0]
    except Exception as e:
        console.print(f"[bold red]Error parsing GenBank file {path}:[/bold red] {str(e)}")
        return None

def parse_all_genbanks(annotations_dir, progress):
    genbank_files = glob.glob(f'{annotations_dir}/*.gbff')
    genbank_paths = {}
    total_files = len(genbank_files)
    task = progress.add_task("[cyan]Parsing genbanks...", total=total_files)
    for i, file in enumerate(genbank_files):
        assembly_contig_id = os.path.basename(file).split('.')[0]
        genbank_paths[assembly_contig_id] = file
        progress.update(task, advance=1, description=f"[cyan]Indexing GenBank files... ({i+1}/{total_files})")
    return genbank_paths

def search_directory(directory):
    alignment_files = []
    for root, dirs, files in os.walk(directory):
        if "align_coords.tsv" in files:
            alignment_files.append(os.path.join(root, "align_coords.tsv"))
    return alignment_files

def find_alignment_files(num_cpus, directory, expected_count, progress):
    start_time = time.time()

    console.print('Getting list of subdirectories!')
    subdirs = [os.path.join(directory, d) for d in os.listdir(directory) if os.path.isdir(os.path.join(directory, d))]
    
    alignment_files = []
    count = 0
    task = progress.add_task("[cyan]Finding alignment files...", total=expected_count)
    with ProcessPoolExecutor(max_workers=num_cpus) as executor:
        console.print('adding jobs to futures!')
        # this is literally unnecessary, we already pull the list of directories above, just tack on the filename. 
        # ugh I hate this why did I write it this way. 
        future_to_dir = {executor.submit(search_directory, subdir): subdir for subdir in subdirs}
        #future_to_dir = {executor.submit(os.path.join, subdir, 'align_coords.tsv'): subdir for subdir in subdirs}
        console.print('All subdirectories added to job queue. Beginning execution.')
        for future in as_completed(future_to_dir):
            subdir = future_to_dir[future]
            try:
                result = future.result()
                alignment_files.extend(result)
                count += 1 
                progress.update(task, advance=1, description=f"[cyan]Finding alignment files... ({count} found)")
                if len(alignment_files) >= expected_count:
                    break
            except Exception as exc:
                console.print(f'Search in {subdir} generated an exception: {exc}')
    console.print('counting alignment files!')
    num_aln_files = len(alignment_files)
    console.print(f'files: {num_aln_files}  count_iter: {count}')
    elapsed_time = time.time() - start_time
    return alignment_files, elapsed_time, num_aln_files
    
def parse_alignment(path):
    with open(path, 'r') as infile:
        lines = infile.readlines()
    keys = lines[0].strip().split('\t')
    return [dict(zip(keys, line.strip().split('\t'))) for line in lines[1:]]

def is_within_range(feature, start, end):
    return int(feature.location.start) >= start and int(feature.location.end) <= end

def get_features_from_range(record, start, end):
    return [feature for feature in record.features if is_within_range(feature, start, end) and feature.type == 'CDS']

def simplify_genes_for_contig(features):
    return [(feature.qualifiers['locus_tag'][0].strip("'[]"), ' '.join(feature.qualifiers['product'])) for feature in features]

def get_coverage(length, coords):
    covered_positions = np.zeros(length, dtype=bool)
    identities = 0
    for start, end, ident in coords:
        if start > end:
            start, end = end, start
        covered_positions[start:end+1] = True
        identities += (end - start + 1) * ident
    percent_coverage = np.sum(covered_positions) / length * 100
    return percent_coverage, identities, covered_positions

def get_homologies(alignments, asm1_record, asm2_record):
    genes_dict = defaultdict(lambda: defaultdict(list))
    for alignment in alignments:
        asm1_name = alignment['QUERY_NAME']
        asm1_start = int(alignment['QUERY_START'])
        asm1_end = int(alignment['QUERY_END'])
        asm1_length = len(asm1_record.seq)
        asm2_name = alignment['REF_NAME']
        asm2_start = int(alignment['REF_START'])
        asm2_end = int(alignment['REF_END'])
        asm2_length = len(asm2_record.seq)
        asm2_identity = float(alignment['IDENTITY'])
        
        asm1_features = get_features_from_range(asm1_record, asm1_start, asm1_end)
        asm2_features = get_features_from_range(asm2_record, asm2_start, asm2_end)
        asm1_genes = simplify_genes_for_contig(asm1_features)
        asm2_genes = simplify_genes_for_contig(asm2_features)
        
        if 'contig_length' not in genes_dict[asm1_name]:
            genes_dict[asm1_name]['contig_length'] = asm1_length
            genes_dict[asm1_name]['assembly_id'] = alignment['QUERY_ID']
        
        genes_dict[asm1_name][asm2_name].append({
            'asm_aln': {
                'asm_name': asm1_name,
                'start': asm1_start,
                'end': asm1_end,
                'aln_length': int(alignment['QUERY_LENGTH']),
                'features': asm1_genes,
            },
            'ref_aln': {
                'ref_id': alignment['REF_ID'],
                'ref_name': asm2_name,
                'start': asm2_start,
                'end': asm2_end,
                'aln_length': int(alignment['REF_LENGTH']),
                'aln_identity': asm2_identity,
                'features': asm2_genes,
                'ref_length': asm2_length,
            },
        })
    return genes_dict

def make_table_for_asm(genes_dict):
    lines = []
    for asm1_id, asm1_data in genes_dict.items():
        assembly1_id = asm1_data['assembly_id']
        asm1_len = asm1_data['contig_length']
        for asm2_id, alignments in asm1_data.items():
            if asm2_id not in ['contig_length', 'assembly_id']:
                asm2_coords = [(cov['asm_aln']['start'], cov['asm_aln']['end'], cov['ref_aln']['aln_identity']) for cov in alignments]
                asm1_aln_coords = [(cov['ref_aln']['start'], cov['ref_aln']['end'], cov['ref_aln']['aln_identity']) for cov in alignments]
                asm2_len = alignments[0]['ref_aln']['ref_length']
                assembly2_id = alignments[0]['ref_aln']['ref_id']
                asm1_genes = [gene for item in alignments for gene in item['ref_aln']['features']]
                asm2_genes = [gene for item in alignments for gene in item['asm_aln']['features']]
                asm1_cov, asm1_identities, _ = get_coverage(asm1_len, asm1_aln_coords)
                asm2_cov, asm2_identities, _ = get_coverage(asm2_len, asm2_coords)
                lines.append([
                    assembly1_id, asm1_id, asm1_len,
                    assembly2_id, asm2_id, asm2_len,
                    f'{asm1_cov:.2f}', asm1_identities, 
                    f'{asm2_cov:.2f}', asm2_identities,
                    len(asm1_genes), len(asm2_genes),
                ])
    return lines

def get_rows(alignments, genes_dict):
    rows = []
    for contig, contig_data in genes_dict.items():
        contig_len = contig_data['contig_length']
        assembly_id = contig_data['assembly_id']
        total_contig_cov = 0
        contig_list = []
        for aln, alignments in contig_data.items():
            if aln not in ['contig_length', 'assembly_id']:
                asm_coords = [(cov['asm_aln']['start'], cov['asm_aln']['end'], cov['ref_aln']['aln_identity']) for cov in alignments]
                ref_coords = [(cov['ref_aln']['start'], cov['ref_aln']['end'], cov['ref_aln']['aln_identity']) for cov in alignments]
                ref_len = alignments[0]['ref_aln']['ref_length']
                contig_cov, _, _ = get_coverage(contig_len, asm_coords)
                ref_cov, _, _ = get_coverage(ref_len, ref_coords)
                contig_list.append(f'({aln} : ({contig_cov:.2f})')
                total_contig_cov += contig_cov
        rows.append([assembly_id, contig, contig_len, f'{total_contig_cov:.2f}', contig_list])
    return rows

def parse_single_pair_aln(alignment_file, genbank_paths):
    asm1_id, asm2_id = os.path.basename(os.path.dirname(alignment_file)).split('_vs_')
    
    if not os.path.getsize(alignment_file) > 1:
        return (False, f'NO HOMOLOGY BETWEEN {asm1_id} and {asm2_id}!', [asm1_id, asm2_id])
    
    asm1_record = parse_genbank(genbank_paths[asm1_id])
    asm2_record = parse_genbank(genbank_paths[asm2_id])
    
    alignments = parse_alignment(alignment_file)
    homologies = get_homologies(alignments, asm1_record, asm2_record)
    detailed_rows = make_table_for_asm(homologies)
    simple_rows = get_rows(alignments, homologies)
    
    return (True, f'Parsed homology between {asm1_id} and {asm2_id}!', simple_rows, detailed_rows)


def parallel_parse(cpus, alignment_files, genbank_paths, progress, expected_alignments):
    simple_rows = [['contig_id','contig_len','total_contig_coverage','list_of_alignments(asm2:percent_cov_of_contig)']]
    detailed_rows = [['assembly1_id', 'asm1_id', 'asm1_len', 'assembly2_id', 'asm2_id', 'asm2_len', 
                      'asm1_cov', 'asm1_identities', 'asm2_cov', 'asm2_identities', 'asm1_genes', 'asm2_genes']]
    no_homology = [['ASM1', 'ASM2']]

    task = progress.add_task("[cyan]Parsing alignments...", total=expected_alignments)
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        console.print('adding jobs to futures!')
        futures = [executor.submit(parse_single_pair_aln, alignment, genbank_paths) for alignment in alignment_files]
        console.print('All jobs added to queue! Beginning parallel execution. Please stand by!')
        for future in as_completed(futures):
            try:
                result = future.result()
                if result[0]:
                    simple_rows.extend(result[2])
                    detailed_rows.extend(result[3])
                else:
                    no_homology.append(result[2])
            except Exception as e:
                console.print(f"[bold red]Error:[/bold red] {e}")
            progress.update(task, advance=1)
    return simple_rows, detailed_rows, no_homology

def generate_alignment_matrix(simple_rows):
    contigs = set()
    for row in simple_rows[1:]:  # Skip header
        contigs.add(row[1])  # contig_id is now in the second column
        for alignment in row[4]:
            contigs.add(alignment.split(':')[0].strip('('))
    
    contigs = sorted(list(contigs))
    matrix = np.zeros((len(contigs), len(contigs)), dtype=float)
    
    for row in simple_rows[1:]:
        source_contig = row[1]  # contig_id is now in the second column
        source_index = contigs.index(source_contig)
        for alignment in row[4]:
            target_contig, coverage = alignment.split(':')
            target_contig = target_contig.strip('(')
            coverage = float(coverage.strip('() '))
            target_index = contigs.index(target_contig)
            matrix[source_index, target_index] = coverage
    
    return matrix, contigs


def write_output_files(output_dir, simple_rows, detailed_rows, no_homologies, progress):
    task = progress.add_task("[cyan]Writing output files...", total=3)
    with open(f'{args.output_dir}/ava_homo_simple.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for row in simple_rows:
            writer.writerow(row)
            progress.update(task, advance=1)
            
    with open(f'{args.output_dir}/ava_homo_detailed.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for row in detailed_rows:
            writer.writerow(row)
            progress.update(task, advance=1)
            
    with open(f'{args.output_dir}/ava_no_homo.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for row in no_homologies:
            writer.writerow(row)
            progress.update(task, advance=1)

def get_args_str(args_dict):
    lines = ''
    for k, v in args_dict.items():
        lines += f'[bold green]{k}[/bold green]: {v}\n'
    return lines

def main():
    parser = argparse.ArgumentParser(description="Homology analysis script")
    parser.add_argument('alignments_dir', help='Directory containing alignments')
    parser.add_argument('annotations_dir', help='Directory containing annotations (genbanks)')
    parser.add_argument('output_dir', help='Directory for outputs')
    parser.add_argument('--cpus', help='Number of CPUs to use', type=int, default=os.cpu_count())
    args = parser.parse_args()
    
    args_dict = vars(args)
    args_str = get_args_str(args_dict)
    
    os.makedirs(args.output_dir, exist_ok=True)
    

            
    progress = Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeRemainingColumn(),
        TimeElapsedColumn(),
    )
    layout = Layout()
    layout.split(
        Layout(Panel(console, title="Console Output", box=box.HEAVY, expand=True), name="top", ratio=2),
        Layout(name="bottom", ratio=1)
    )
    layout["bottom"].split_row(
            Layout(Panel(progress, title="Task Progress", box=box.HEAVY), name="bottom_right", ratio=1),
            Layout(Panel(args_str, title="Arguments", box=box.HEAVY,), name="bottom_left", ratio=1),
    )
        
    with Live(layout, refresh_per_second=10) as live:
        console.print(Panel.fit("[bold cyan]Starting Homology Analysis[/bold cyan]", highlight=True, box=box.ROUNDED))

        console.print("[yellow]Indexing all GenBank files...[/yellow]")
        genbank_paths = parse_all_genbanks(args.annotations_dir, progress)
        
        if not genbank_paths:
            console.print("[bold red]Error: No valid GenBank files found. Exiting.[/bold red]")
            sys.exit(1)

        n_genbanks = len(genbank_paths)
        console.print(f"[green]Successfully indexed {n_genbanks} GenBank files.[/green]")
        console.print()  # Add an empty line for spacing

        expected_alignments = (n_genbanks * (n_genbanks - 1)) // 2
        console.print(f"[yellow]Expected number of alignment files: {expected_alignments}[/yellow]")

        console.print("[yellow]Finding alignment files...[/yellow]")
        alignment_files, elapsed_time, num_aln_files = find_alignment_files(args.cpus, args.alignments_dir, expected_alignments, progress)

        console.print(f"[green]Found {len(alignment_files)} alignment files in {elapsed_time:.2f} seconds.[/green]")

        if num_aln_files < expected_alignments:
            console.print(f"[bold yellow]Warning: Found fewer alignment files than expected. Expected: {expected_alignments}, Found: {len(alignment_files)}[/bold yellow]")
        console.print()
        console.print("[yellow]Starting parallel parsing...[/yellow]")
        simple_rows, detailed_rows, no_homologies = parallel_parse(args.cpus, alignment_files, genbank_paths, progress, num_aln_files)

        console.print()
        console.print("[yellow]Writing output files...[/yellow]")
        #writing_task = progress.add_task("[cyan]Writing output files...", total=3)
        write_output_files(args.output_dir, simple_rows, detailed_rows, no_homologies, progress)#, writing_task)
        #progress.update(writing_task, completed=3)

        console.print()
        console.print("[yellow]Generating alignment matrix...[/yellow]")
        matrix_task = progress.add_task("[cyan]Generating alignment matrix...", total=1)
        matrix, contigs = generate_alignment_matrix(simple_rows)
        progress.update(matrix_task, completed=1)

        console.print()
        console.print("[yellow]Saving alignment matrix...[/yellow]")
        saving_task = progress.add_task("[cyan]Saving alignment matrix...", total=1)
        np.savetxt(f'{args.output_dir}/alignment_matrix.csv', matrix, delimiter=',', header=','.join(contigs), comments='')
        progress.update(saving_task, completed=1)
        
        console.print()
        console.print(Panel.fit("[bold green]Homology Analysis Completed Successfully![/bold green]"))

if __name__ == "__main__":
    main()