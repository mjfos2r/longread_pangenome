import os
import sys
import csv
import glob
import re

from Bio import SeqIO
from collections import defaultdict
import argparse
from pathlib import Path
import traceback
import pickle

import pandas as pd
import numpy as np

import logging
from logging.handlers import BufferingHandler

import tempfile
import shutil
import time
from datetime import datetime

from multiprocessing import cpu_count
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed, wait, FIRST_COMPLETED
from queue import Queue, Empty
from threading import Thread, Lock

from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeRemainingColumn, TimeElapsedColumn, SpinnerColumn
from rich.console import Console, RenderableType
from rich.panel import Panel
from rich.live import Live
from rich.table import Table
from rich.layout import Layout
from rich.text import Text
from rich import box
from rich.align import Align

from collections import deque

# This is to get correct log printing behavior in our pane!
# instead of flushing the entire list once the capacity has been reached, it simply pops the oldest msg from the list.
class RollingBufferHandler(logging.Handler):
    def __init__(self, capacity):
        super().__init__()
        self.capacity = capacity
        self.buffer = deque(maxlen=capacity)
        self.queue = Queue()
        self.console = Console(record=True)

    def emit(self, record):
        try:
            msg = self.format(record)
            rich_text = self.console.render_str(msg)
            self.queue.put(rich_text)
            self._update_buffer()
        except Exception:
            self.handleError(record)

    def _update_buffer(self):
        while not self.queue.empty():
            try:
                msg = self.queue.get_nowait()
                if len(self.buffer) == self.capacity:
                    self.buffer.popleft()
                self.buffer.append(msg)
            except Empty:
                break

    def get_logs(self):
        self._update_buffer()
        return self.buffer

class ConsolePanel(Console):
    def __init__(self,*args,**kwargs):
        console_file = open(os.devnull,'w')
        super().__init__(record=True,file=console_file,*args,**kwargs)

    def __rich_console__(self,console,options,):#stderr=False):
        texts = self.export_text(clear=False,).split('\n') # styles=True <- This fucks up the borders.
        for line in texts[-options.height:]:
            yield line

class LogPanel:
    def __init__(self, RollingBufferHandler):
        self.buffer_handler = RollingBufferHandler

    def __rich__(self) -> Panel:
        log_content = Text('\n').join(self.buffer_handler.get_logs())
        return Panel(log_content, box=box.ROUNDED, title="Recent Messages", expand=True)

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
    task = progress.add_task("[cyan]Parsing genbanks...", total=total_files, expand=True)
    for i, file in enumerate(genbank_files):
        assembly_contig_id = os.path.basename(file).split('.')[0]
        genbank_paths[assembly_contig_id] = file
        progress.update(task, advance=1, description=f"[cyan]Indexing GenBank files... ({i+1}/{total_files})", expand=True)
    return genbank_paths

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
        identities += (end - start + 1) * float(ident/100)
    total_homology = np.sum(covered_positions)
    percent_coverage = total_homology / length * 100
    return percent_coverage, identities, total_homology

class CursedDict(defaultdict):
    def __init__(self):
        super().__init__(list)
        self["assembly_id"] = ""
        self["contig_length"] = int()

    def __getitem__(self, key):
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        if key == "assembly_id" and not isinstance(value, str):
            raise ValueError("assembly_id must be a string")
        elif key == "contig_length" and not isinstance(value, int):
            raise ValueError("contig_length must be an integer")
        super().__setitem__(key, value)

def get_homologies(alignments, genbank_paths):
    genes_dict = defaultdict(CursedDict)
    for alignment in alignments:
        if isinstance(alignment, str):
            # Single file processing
            alignment_data = alignment.strip().split('\n')
            keys = alignment_data[0].strip().split('\t')
            values = alignment_data[1].strip().split('\t')
            alignment_dict = dict(zip(keys, values))
        else:
            # Directory processing
            alignment_dict = alignment

        asm1_id = alignment_dict['QUERY_ID']
        asm1_name = alignment_dict['QUERY_NAME']
        asm1_start = int(alignment_dict['QUERY_START'])
        asm1_end = int(alignment_dict['QUERY_END'])
        asm1_length = int(len(parse_genbank(genbank_paths[asm1_id]).seq))
        asm2_id = alignment_dict['REF_ID']
        asm2_name = alignment_dict['REF_NAME']
        asm2_start = int(alignment_dict['REF_START'])
        asm2_end = int(alignment_dict['REF_END'])
        asm2_length = int(len(parse_genbank(genbank_paths[asm2_id]).seq))
        asm2_identity = float(alignment_dict['IDENTITY'])

        asm1_record = parse_genbank(genbank_paths[asm1_id])
        asm2_record = parse_genbank(genbank_paths[asm2_id])
        asm1_features = get_features_from_range(asm1_record, asm1_start, asm1_end)
        asm2_features = get_features_from_range(asm2_record, asm2_start, asm2_end)
        asm1_genes = simplify_genes_for_contig(asm1_features)
        asm2_genes = simplify_genes_for_contig(asm2_features)

        genes_dict[asm1_id]['contig_length'] = asm1_length
        genes_dict[asm1_id]['assembly_id'] = alignment_dict['QUERY_ID']

        if asm2_id not in genes_dict[asm1_id]:
            genes_dict[asm1_id][asm2_id] = []

        genes_dict[asm1_id][asm2_id].append({
            'asm_aln': {
                'asm_name': asm1_name,
                'start': asm1_start,
                'end': asm1_end,
                'aln_length': int(alignment_dict['QUERY_LENGTH']),
                'features': asm1_genes,
            },
            'ref_aln': {
                'ref_id': asm2_id,
                'ref_name': asm2_name,
                'start': asm2_start,
                'end': asm2_end,
                'aln_length': int(alignment_dict['REF_LENGTH']),
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
                asm1_aln_coords = [(cov['asm_aln']['start'], cov['asm_aln']['end'], cov['ref_aln']['aln_identity']) for cov in alignments]
                asm2_aln_coords = [(cov['ref_aln']['start'], cov['ref_aln']['end'], cov['ref_aln']['aln_identity']) for cov in alignments]
                asm2_len = alignments[0]['ref_aln']['ref_length']
                assembly2_id = alignments[0]['ref_aln']['ref_id']
                asm2_genes = [gene for item in alignments for gene in item['ref_aln']['features']]
                asm1_genes = [gene for item in alignments for gene in item['asm_aln']['features']]
                asm1_cov, asm1_identities, asm1_total_homology = get_coverage(asm1_len, asm1_aln_coords)
                asm2_cov, asm2_identities, asm2_total_homology = get_coverage(asm2_len, asm2_aln_coords)
                lines.append([
                    assembly1_id, asm1_id, asm1_len,
                    assembly2_id, asm2_id, asm2_len,
                    asm1_cov, asm1_total_homology, asm1_identities,
                    asm2_cov, asm2_total_homology, asm2_identities,
                    len(asm1_genes), len(asm2_genes),
                ])
    return lines

def get_rows(alignments, genes_dict):
    rows = []
    for contig, contig_data in genes_dict.items():
        contig_len = contig_data['contig_length']
        assembly_id = contig_data['assembly_id']
        total_contig_cov = 0
        contig_alignments = {}

        if isinstance(contig_len, list):
            contig_len = contig_len[0]

        for aln, alignments in contig_data.items():
            if aln not in ['contig_length', 'assembly_id']:
                asm_coords = [(cov['asm_aln']['start'], cov['asm_aln']['end'], cov['ref_aln']['aln_identity']) for cov in alignments]
                ref_coords = [(cov['ref_aln']['start'], cov['ref_aln']['end'], cov['ref_aln']['aln_identity']) for cov in alignments]
                ref_len = alignments[0]['ref_aln']['ref_length']
                _, _, contig_total_homology = get_coverage(contig_len, asm_coords)
                _, _, ref_total_homology = get_coverage(ref_len, ref_coords)
                contig_alignments[aln] = int(contig_total_homology)
                total_contig_cov += contig_total_homology

        if contig_alignments and contig_len > 0:
            rows.append([
                assembly_id,
                contig,
                int(contig_len),
                int(total_contig_cov),
                [f"({k}:{v})" for k, v in contig_alignments.items()]
            ])
    return rows

def check_tsv_content(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        return True if len(lines) > 1 else False

def process_alignment_nucl(alignment, genbank_paths):
    if isinstance(alignment, str):
        # Single file processing
        alignment_data = alignment.strip().split('\n')
        keys = alignment_data[0].strip().split('\t')
        values = alignment_data[1].strip().split('\t')
        alignment_dict = dict(zip(keys, values))
    else:
        # Directory processing
        alignment_dict = alignment

    asm1_id = alignment_dict['QUERY_ID']
    asm2_id = alignment_dict['REF_ID']

    asm1_record = parse_genbank(genbank_paths[asm1_id])
    asm2_record = parse_genbank(genbank_paths[asm2_id])

    if not asm1_record or not asm2_record:
        simple_row = [asm1_id, asm1_id, len(asm1_record.seq) if asm1_record else 0, '0', []]
        detailed_row = [asm1_id, asm1_id, len(asm1_record.seq) if asm1_record else 0,
                        asm2_id, asm2_id, len(asm2_record.seq) if asm2_record else 0,
                        0, 0, 0, 0, 0, 0, 0, 0]
        return (True, f'NO HOMOLOGY BETWEEN {asm1_id} and {asm2_id}!', [simple_row], [detailed_row])

    homologies = get_homologies([alignment_dict], genbank_paths)
    detailed_rows = make_table_for_asm(homologies)
    simple_rows = get_rows([alignment_dict], homologies)

    return (True, f'Parsed homology between {asm1_id} and {asm2_id}!', simple_rows, detailed_rows)

def process_alignment_gene(alignment, genbank_paths):
    if isinstance(alignment, str):
        # Single file processing
        alignment_data = alignment.strip().split('\n')
        keys = alignment_data[0].strip().split('\t')
        values = alignment_data[1].strip().split('\t')
        alignment_dict = dict(zip(keys, values))
    else:
        # Directory processing
        alignment_dict = alignment

    gene_count = defaultdict()
    asm1_id = alignment_dict['QUERY_ID']
    asm2_id = alignment_dict['REF_ID']

    asm1_record = parse_genbank(genbank_paths[asm1_id])
    asm2_record = parse_genbank(genbank_paths[asm2_id])
    
    gene_count[asm1_id] = len([feature for feature in asm1_record.features 
                               if feature.type == "CDS"])
    gene_count[asm2_id] = len([feature for feature in asm2_record.features 
                               if feature.type == "CDS"])

    gene_homology_rows = get_gene_homology([alignment_dict], genbank_paths)

    return (True, f'Parsed gene homology between {asm1_id} and {asm2_id}!', gene_homology_rows, None, gene_count)

def get_gene_homology(alignments, genbank_paths):
    gene_homology_rows = []
    for alignment in alignments:
        if isinstance(alignment, str):
            # Single file processing
            alignment_data = alignment.strip().split('\n')
            keys = alignment_data[0].strip().split('\t')
            values = alignment_data[1].strip().split('\t')
            alignment_dict = dict(zip(keys, values))
        else:
            # Directory processing
            alignment_dict = alignment

        asm1_id = alignment_dict['QUERY_ID']
        asm1_start = int(alignment_dict['QUERY_START'])
        asm1_end = int(alignment_dict['QUERY_END'])
        asm2_id = alignment_dict['REF_ID']
        asm2_start = int(alignment_dict['REF_START'])
        asm2_end = int(alignment_dict['REF_END'])
        identity = float(alignment_dict['IDENTITY'])

        asm1_record = parse_genbank(genbank_paths[asm1_id])
        asm2_record = parse_genbank(genbank_paths[asm2_id])

        asm1_features = get_features_from_range(asm1_record, asm1_start, asm1_end)
        asm2_features = get_features_from_range(asm2_record, asm2_start, asm2_end)

        for asm1_feature, asm2_feature in zip(asm1_features, asm2_features):
            gene_homology_row = [
                asm1_id,
                asm1_feature.qualifiers['locus_tag'][0],
                asm1_feature.qualifiers.get('product', [''])[0],
                f"{asm1_start}-{asm1_end}",
                asm2_id,
                asm2_feature.qualifiers['locus_tag'][0],
                asm2_feature.qualifiers.get('product', [''])[0],
                f"{asm2_start}-{asm2_end}",
                identity
            ]
            gene_homology_rows.append(gene_homology_row)

    return gene_homology_rows

def discover_files(directory, queue):
    discovered_count = 0
    logging.debug(f'discover_files: Beginning walk through {directory}')
    for dirpath, _, filenames in os.walk(directory):
        for filename in filenames:
            if filename == "align_coords.tsv":
                full_path = os.path.join(dirpath, filename)
                queue.put(full_path)
                discovered_count += 1
    logging.debug('discover_files: inserting None into queue')
    queue.put(None)
    logging.info(f"[bold green]Discovery complete. Found[/bold green] {discovered_count} [bold green]files[/bold green]")

def worker(file_list, genbank_paths, analysis_type):
    results = []
    with ThreadPoolExecutor() as thread_executor:
        if analysis_type == 'nucl':
            futures = [thread_executor.submit(process_alignment_nucl, file_path, genbank_paths) for file_path in file_list]
        else:  # gene homology
            futures = [thread_executor.submit(process_alignment_gene, file_path, genbank_paths) for file_path in file_list]
        for future in as_completed(futures):
            results.append(future.result())
    return results

def find_and_process_alignments(directory_or_file, genbank_paths, num_cpus, progress, expected_count, max_queue_size, analysis_type, all_alns_file=None):
    start_time = time.time()

    if analysis_type == 'nucl':
        simple_rows = [['assembly_id', 'contig_id', 'contig_len', 'total_contig_homology', 'list_of_alignments(asm2:homology_length)']]
        detailed_rows = [['assembly1_id', 'asm1_id', 'asm1_len', 'assembly2_id', 'asm2_id', 'asm2_len',
                          'asm1_cov', 'asm1_total_homology', 'asm1_identities', 'asm2_cov', 'asm2_total_homology', 'asm2_identities', 'asm1_genes', 'asm2_genes']]
    else:  # gene homology
        gene_homology_rows = [['asm1_id', 'asm1_locus_tag', 'asm1_product', 'asm1_coords', 'asm2_id', 'asm2_locus_tag', 'asm2_product', 'asm2_coords', 'identity']]

    empty_alignment_files = []
    gene_count = defaultdict(int)

    logging.debug('[bold blue]Beginning alignment processing[/bold blue]')

    if all_alns_file:
        # Process single file
        logging.info(f"Processing single alignment file: {all_alns_file}")
        processing_task = progress.add_task("[cyan]Processing alignment file...", total=None, expand=True)
        
        with open(all_alns_file, 'r') as f:
            alignments = f.read().split('\n\n')  # Assuming alignments are separated by blank lines
        
        num_alignments = len(alignments)
        progress.update(processing_task, total=num_alignments)
        
        # Create a structure similar to what directory-based processing would use
        alignment_dict = {}
        for alignment in alignments:
            alignment_data = alignment.strip().split('\n')
            keys = alignment_data[0].strip().split('\t')
            values = alignment_data[1].strip().split('\t')
            aln_dict = dict(zip(keys, values))
            asm1_id = aln_dict['QUERY_ID']
            asm2_id = aln_dict['REF_ID']
            if asm1_id not in alignment_dict:
                alignment_dict[asm1_id] = {}
            alignment_dict[asm1_id][asm2_id] = aln_dict

        # Process alignments
        processed_count = 0
        for asm1_id, asm2_dict in alignment_dict.items():
            for asm2_id, aln_dict in asm2_dict.items():
                if analysis_type == 'nucl':
                    result = process_alignment_nucl(aln_dict, genbank_paths)
                    simple_rows.extend(result[2])
                    detailed_rows.extend(result[3])
                else:  # gene homology
                    result = process_alignment_gene(aln_dict, genbank_paths)
                    gene_homology_rows.extend(result[2])
                    gene_count.update(result[4])
                processed_count += 1
                if processed_count % 1000 == 0:
                    progress.update(processing_task, advance=1000, description=f"[cyan]Processed {processed_count}/{num_alignments} alignments")

        progress.update(processing_task, advance=num_alignments % 1000, description=f"[cyan]Processed all {num_alignments} alignments")

    else:
        # Existing directory-based processing logic
        logging.debug('[yellow]Beginning file discovery, initializing file Queue() and discovery Thread()[/yellow]')
        discovery_task = progress.add_task("[cyan]Discovering alignment files...", total=None, expand=True)
        processing_task = progress.add_task("[cyan]Processing alignment files...", total=expected_count, expand=True)
        
        file_queue = Queue()
        discovery_thread = Thread(target=discover_files, args=(directory_or_file, file_queue), name="DiscoveryThread")

        logging.debug('Starting discovery thread!')
        discovery_thread.start()

        discovered_count = 0
        processed_count = 0
        
        with ProcessPoolExecutor(max_workers=num_cpus) as process_executor:
            futures = []
            discovery_complete = False
            all_files_processed = False
            file_buffer = []

            while not discovery_complete or futures or file_buffer:
                while not discovery_complete and len(file_buffer) < num_cpus * 1000:
                    try:
                        aln_path = file_queue.get(timeout=1)
                        if aln_path is None:
                            discovery_complete = True
                            logging.debug("`None` path reached! Discovery is finished. Waiting for discovery thread to complete. Merging thread with main!")
                            discovery_thread.join()
                            logging.debug("Discovery thread successfully merged!")
                            progress.update(discovery_task, completed=discovered_count, total=discovered_count,
                                            description=f"[bold green]File discovery complete! Found[/bold green] {discovered_count} [bold green]alignments.[/bold green]"
                                            )
                            break
                        else:
                            file_buffer.append(aln_path)
                            discovered_count += 1
                            progress.update(discovery_task, advance=discovered_count,
                                            description=f"[cyan]Discovering alignment files... ({discovered_count} found!)"
                                            )
                    except Empty:
                        if not discovery_thread.is_alive():
                            discovery_complete = True
                            logging.debug("discovery thread is dead and queue is empty. merging the corpse with main thread.")
                            discovery_thread.join()
                            logging.debug("discovery thread merged with main.")
                            progress.update(discovery_task, completed=discovered_count, total=discovered_count, description=f"[bold green]File discovery complete! Found[/bold green] {discovered_count} [bold green]alignments.[/bold green]")
                            break
                logging.debug("Chunksize reached. Spawning jobs!")
                while file_buffer and len(futures) < num_cpus:
                    chunk_size = len(file_buffer) // num_cpus + 1
                    file_chunk = file_buffer[:chunk_size]
                    file_buffer = file_buffer[chunk_size:]
                    futures.append(process_executor.submit(worker, file_chunk, genbank_paths, analysis_type))

                for future in as_completed(futures):
                    futures.remove(future)
                    try:
                        results = future.result()
                        for result in results:
                            if not result[0]:  # If alignment file is empty
                                empty_alignment_files.append(result[-1])
                            if analysis_type == 'nucl':
                                simple_rows.extend(result[2])
                                detailed_rows.extend(result[3])
                            else:  # gene homology
                                gene_homology_rows.extend(result[2])
                                gene_count.update(result[4])
                        processed_count += len(results)
                        progress.update(
                            processing_task, advance=len(results),
                            description=f"[cyan]Processing alignment files... ({processed_count} processed)"
                        )
                    except Exception as e:
                        logging.error(f"[bold red]Error processing file:[/bold red]: {traceback.format_exc()}")

                    if processed_count >= expected_count and not all_files_processed:
                        logging.info(
                            f"Reached expected count of {expected_count}, setting [white]all_files_processed[/white] flag to [green]True[/green]!")
                        all_files_processed = True
                        progress.update(processing_task, total=discovered_count)

    elapsed_time = time.time() - start_time
    logging.info(f'[bold green]Processing complete![/bold green]\nTotal Runtime: {elapsed_time:.2f} seconds. [bold green]Alignments processed:[/bold green] {processed_count}')

    if analysis_type == 'nucl':
        logging.info(f"Total rows in simple_rows: {len(simple_rows) - 1}")
        df = pd.DataFrame(simple_rows[1:], columns=simple_rows[0])
        logging.info(f"Total rows in DataFrame: {len(df)}")
        logging.debug('aggregating alignments by contig')
        grouped_data = aggregate_alignments(df)
        logging.info(f"Total contigs in grouped_data: {len(grouped_data)}")
        logging.debug('creating summary dataframe of aggregated alignments')
        summary_df = create_summary_dataframe(grouped_data)
        logging.info(f"Total rows in summary DataFrame: {len(summary_df)}")
        logging.debug('creating matrices from aggregated alignments!')
        np_matrix, df_matrix, contig_labels = create_homology_matrix(grouped_data)
        return simple_rows, detailed_rows, summary_df, np_matrix, df_matrix, contig_labels, processed_count, elapsed_time, empty_alignment_files
    else:  # gene homology
        return gene_homology_rows, processed_count, elapsed_time, empty_alignment_files, gene_count

def aggregate_alignments(df):
    grouped = {}
    all_contigs = set()

    for _, row in df.iterrows():
        source_contig = row['contig_id']
        all_contigs.add(source_contig)

        if source_contig not in grouped:
            grouped[source_contig] = {
                'source_len': row['contig_len'],
                'total_contig_homology': row['total_contig_homology'],
                'alignments': {}
            }

        alignments = row['list_of_alignments(asm2:homology_length)']
        if isinstance(alignments, str):
            alignments = eval(alignments)  # Be cautious with eval, ensure input is safe

        for alignment in alignments:
            match = re.match(r"\(([^:]+):(\d+)\)", alignment)
            if match:
                target_contig, homology_len = match.groups()
                homology_len = int(homology_len)
                all_contigs.add(target_contig)
                grouped[source_contig]['alignments'][target_contig] = homology_len

                if target_contig not in grouped:
                    grouped[target_contig] = {
                        'source_len': None,
                        'total_contig_homology': None,
                        'alignments': {}
                    }

    for contig in all_contigs:
        if contig not in grouped:
            grouped[contig] = {
                'source_len': None,
                'total_contig_homology': None,
                'alignments': {}
            }

    return grouped

def create_summary_dataframe(grouped_data):
    summary_data = []
    logging.debug('beginning summary df construction!')
    for source_contig, data in grouped_data.items():
        summary_data.append({
            'source_contig': source_contig,
            'source_len': data['source_len'],
            'total_contig_homology': data['total_contig_homology'],
            'num_alignments': len(data['alignments']),
            'target_contigs': ', '.join(data['alignments'].keys()),
            'homology_lengths': ', '.join(map(str, data['alignments'].values()))
        })
    logging.debug('finished summary df construction!')
    return pd.DataFrame(summary_data)

def create_homology_matrix(grouped_data):
    all_contigs = list(grouped_data.keys())
    n = len(all_contigs)
    contig_to_index = {contig: i for i, contig in enumerate(all_contigs)}
    np_matrix = np.zeros((n, n), dtype=np.int64)

    for source_contig, data in grouped_data.items():
        source_index = contig_to_index[source_contig]
        source_len = data.get('source_len')
        if source_len is not None:
            np_matrix[source_index, source_index] = source_len

        alignments = data.get('alignments', {})
        for target_contig, homology_len in alignments.items():
            target_index = contig_to_index[target_contig]
            np_matrix[source_index, target_index] = homology_len
            np_matrix[target_index, source_index] = homology_len  # Symmetric
    df_matrix = pd.DataFrame(np_matrix, index=all_contigs, columns=all_contigs)
    return np_matrix, df_matrix, all_contigs

def create_shared_gene_matrix(gene_homology_rows, identity_cutoff, gene_count):
    # Extract assembly IDs
    asm_ids = sorted(set(row[0] for row in gene_homology_rows[1:]) | set(row[4] for row in gene_homology_rows[1:]))
    n_asm = len(asm_ids)

    # Create a dictionary to map assembly IDs to matrix indices
    asm_to_index = {asm: i for i, asm in enumerate(asm_ids)}

    # Initialize the matrix with zeros
    shared_gene_matrix = np.zeros((n_asm, n_asm), dtype=int)

    # Populate the matrix
    for row in gene_homology_rows[1:]:  # Skip header row
        if float(row[-1]) >= identity_cutoff:
            asm1_idx = asm_to_index[row[0]]
            asm2_idx = asm_to_index[row[4]]
            shared_gene_matrix[asm1_idx, asm2_idx] += 1
            shared_gene_matrix[asm2_idx, asm1_idx] += 1  # Symmetric matrix
            shared_gene_matrix[asm1_idx, asm1_idx] = gene_count[asm1_idx]
            shared_gene_matrix[asm2_idx, asm2_idx] = gene_count[asm2_idx] # add the gene count to diagonal
            
    return pd.DataFrame(shared_gene_matrix, index=asm_ids, columns=asm_ids)

def write_output_files_nucl(output_dir, simple_rows, detailed_rows, summary_df, np_matrix, df_matrix, contig_labels, progress, version, empty_alignment_files):
    logging.info('[bold_blue]Writing rows to files. Please stand by.[/bold_blue]')
    task = progress.add_task("[cyan]Writing output files...", total=7)

    with open(f'{output_dir}/ava_homology_simple_{version}.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for row in simple_rows:
            writer.writerow(row)
    progress.update(task, advance=1)
    logging.info(f'[bold_green]simple_rows file has been written![/bold_green]')

    with open(f'{output_dir}/ava_homology_detailed_{version}.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for row in detailed_rows:
            writer.writerow(row)
    progress.update(task, advance=1)
    logging.info(f'[bold_green]detailed_rows file has been written![/bold_green]')

    summary_df.to_csv(f'{output_dir}/ava_homology_summary_{version}.tsv', sep='\t', index=False)
    progress.update(task, advance=1)
    logging.info(f'[bold_green]summary DataFrame has been written![/bold_green]')

    np.savetxt(f'{output_dir}/ava_homology_matrix_{version}.csv', np_matrix, delimiter=',')
    progress.update(task, advance=1)
    logging.info(f'[bold_green]numpy matrix has been written![/bold_green]')

    df_matrix.to_csv(f'{output_dir}/ava_homology_matrix_labeled_{version}.csv')
    progress.update(task, advance=1)
    logging.info(f'[bold_green]pandas DataFrame matrix has been written![/bold_green]')

    with open(f'{output_dir}/ava_homology_np_matrix_labels_{version}.txt', 'w') as f:
        for label in contig_labels:
            f.write(f"{label}\n")
    progress.update(task, advance=1)

        # Write empty alignment files
    with open(f'{output_dir}/empty_alignment_files_{version}.txt', 'w') as f:
        for file in empty_alignment_files:
            f.write(f"{file}\n")
    progress.update(task, advance=1)
    logging.info(f'[bold_green]Empty alignment files list has been written![/bold_green]')

    logging.info(f'[bold_green]contig labels have been written![/bold_green]')

    logging.info(f'[bold_green]All output files have been written![/bold_green]')

def write_output_files_gene(output_dir, gene_homology_rows, progress, version, min_identity, identity_step, empty_alignment_files):
    logging.info('[bold_blue]Writing rows to files. Please stand by.[/bold_blue]')

    num_files = len(range(int(min_identity), 101, int(identity_step))) * 2 + 2  # Added 1 for empty alignment files
    task = progress.add_task("[cyan]Writing output files...", total=num_files)

    logging.debug('writing full gene homology table!')
    with open(f'{output_dir}/gene_homology_full_{version}.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for row in gene_homology_rows:
            writer.writerow(row)
    progress.update(task, advance=1)
    logging.info(f'[bold_green]Full gene homology table has been written![/bold_green]')

    # Separate header and data rows
    header = gene_homology_rows[0]
    data_rows = gene_homology_rows[1:]

    os.makedirs(f'{output_dir}/tables', exist_ok=True)
    os.makedirs(f'{output_dir}/matrices', exist_ok=True)
    for identity_cutoff in range(int(min_identity), 101, int(identity_step)):
        filtered_rows = [row for row in data_rows if float(row[-1]) >= identity_cutoff]

        with open(f'{output_dir}/tables/gene_homology_{identity_cutoff}percent_{version}.tsv', 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(header)  # Write header
            for row in filtered_rows:
                writer.writerow(row)

        progress.update(task, advance=1)
        logging.info(f'[bold_green]Gene homology table for {identity_cutoff}% identity cutoff has been written![/bold_green]')

        # Generate and write shared gene matrix
        shared_gene_matrix = create_shared_gene_matrix([header] + filtered_rows, identity_cutoff)
        shared_gene_matrix.to_csv(f'{output_dir}/matrices/shared_gene_matrix_{identity_cutoff}percent_{version}.tsv', sep='\t')

        progress.update(task, advance=1)
        logging.info(f'[bold_green]Shared gene matrix for {identity_cutoff}% identity cutoff has been written![/bold_green]')

    # Write empty alignment files
    with open(f'{output_dir}/empty_alignment_files_{version}.txt', 'w') as f:
        for file in empty_alignment_files:
            f.write(f"{file}\n")
    progress.update(task, advance=1)
    logging.info(f'[bold_green]Empty alignment files list has been written![/bold_green]')

    logging.info(f'[bold_green]All output files have been written![/bold_green]')

def setup_logging(log_file):
    logger = logging.getLogger()
    logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
    )
    return logger

def create_layout(console, progress, args_str, buffering_handler):
    log_panel = LogPanel(buffering_handler)

    layout = Layout()
    layout.split(
        Layout(name="top", ratio=3),
        Layout(name="bottom", ratio=1)
    )
    layout["top"].split_row(
        Layout(Panel(console, title="Console Output", box=box.HEAVY, expand=True), name="top_left", ratio=1),
        Layout(Panel(log_panel, title="DEBUG", box=box.HEAVY, expand=True), name="top_right", ratio=2),
    )
    layout["bottom"].split_row(
        Layout(Panel(progress, title="Task Progress", box=box.HEAVY,), name="bottom_right", ratio=2),
        Layout(Panel.fit(args_str, title="Arguments", box=box.HEAVY), name="bottom_left", ratio=1),
    )
    return layout

def parse_arguments():
    parser = argparse.ArgumentParser(description="Homology Analysis")
    parser.add_argument('--annotations_dir', required=True, help="Directory containing GenBank files")
    parser.add_argument('--alignments_dir', required=False, help="Directory containing alignment files")
    parser.add_argument('--all_alns_file', required=False, help="Single file containing all concatenated alignments")
    parser.add_argument('--output_dir', required=True, help="Directory to store output files")
    parser.add_argument('--logfile', required=True, help="Path to log file")
    parser.add_argument('--cores', type=int, help="Number of CPU cores to use", required=False)
    parser.add_argument('--dale', action='store_true', help="Redline this server by using all cores (less two for headroom)", required=False)
    parser.add_argument('--version', type=str, help="Which version is this?", required=True)
    parser.add_argument('--analysis_type', choices=['nucl', 'gene'], required=True, help="Type of analysis to perform")
    parser.add_argument('--min_identity', type=float, default=50.0, help="Minimum identity cutoff for gene homology (default: 50.0)")
    parser.add_argument('--identity_step', type=float, default=10.0, help="Identity increment step for gene homology (default: 10.0)")
    
    args = parser.parse_args()
    
    if not args.alignments_dir and not args.all_alns_file:
        parser.error("Either --alignments_dir or --all_alns_file must be provided")
    if args.alignments_dir and args.all_alns_file:
        parser.error("Only one of --alignments_dir or --all_alns_file should be provided")
    
    return args

def main(args):
    console = ConsolePanel()

    temp_log_file = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.log')
    temp_log_path = temp_log_file.name
    temp_log_file.close()
    try:
        main_logger = setup_logging(temp_log_path)
        main_logger.setLevel(logging.DEBUG)
        logging.debug('initialized logger: main_logger')

        LOG_BUFFER_MAX_MSGS = 16
        logging.debug(f"max buffer length: {LOG_BUFFER_MAX_MSGS}")
        logging.debug("setting up log formatter for pane view")
        formatter =  logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        logging.debug("setting up custom rolling buffering handler for pane view")
        buffering_handler = RollingBufferHandler(capacity=LOG_BUFFER_MAX_MSGS)
        buffering_handler.setFormatter(formatter)
        logging.debug("adding custom rolling buffering handler to main_logger")
        main_logger.addHandler(buffering_handler)
        progress = Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TaskProgressColumn(),
            TimeRemainingColumn(),
            TimeElapsedColumn(),
            SpinnerColumn(),
        )
        logging.debug("parsing args to dict and creating string for pane view")
        args_dict = vars(args)
        args_str  = "\n".join(f"[bold green]{k}[/bold green]: {v}" for k,v in args_dict.items())
        logging.debug(f"creating path for subdirectory:{args.version} in output:{args.output_dir}")
        args.output_dir = os.path.join(args.output_dir, args.version)
        logging.debug(f"new output path: {args.output_dir}")

        layout = create_layout(console, progress, args_str, buffering_handler)

        logging.debug("Starting live view")
        with Live(layout, refresh_per_second=10) as live:
            logging.info(f"Using the following arguments for this run\n{args_str}\n")

            logging.info(f'Creating output directory: {args.output_dir}')
            console.print(f'Creating output directory: {args.output_dir}')
            os.makedirs(args.output_dir, exist_ok=True)
            logging.info('[green]Output directory successfully created[/green]')
            console.print('[green]Output directory successfully created[/green]')

            console.print(Panel.fit("[bold cyan]Starting Homology Analysis[/bold cyan]", highlight=True, box=box.ROUNDED))
            logging.info("[bold cyan]Starting Homology Analysis[/bold cyan]")

            console.print("[bold blue]Indexing all GenBank files...[/bold blue]")
            logging.info("[bold blue]Indexing all GenBank files...[/bold blue]")
            genbank_paths = parse_all_genbanks(args.annotations_dir, progress)

            if not genbank_paths:
                console.print("[bold red]Error: No valid GenBank files found. Exiting.[/bold red]")
                logging.error("[bold red]Error: No valid GenBank files found. Exiting.[/bold red]")
                sys.exit(1)

            n_genbanks = len(genbank_paths)
            console.print(f"[green]Successfully indexed[/green] [green]{n_genbanks}[/green] [green]GenBank files.[/green]")
            logging.info(f"[green]Successfully indexed[/green] [green]{n_genbanks}[/green] [green]GenBank files.[/green]")

            expected_alignments = (n_genbanks * (n_genbanks - 1)) // 2
            console.print(f"[yellow]Expected number of alignment files:[/yellow] {expected_alignments}")
            logging.info(f"[yellow]Expected number of alignment files:[/yellow] {expected_alignments}")
            live.refresh()

            console.print("[bold blue]Finding and processing alignment files...[/bold blue]")
            logging.info("[bold blue]Finding and processing alignment files...[/bold blue]")
            max_queue_size = 10000

            if args.all_alns_file:
                if args.analysis_type == 'nucl':
                    simple_rows, detailed_rows, summary_df, np_matrix, df_matrix, contig_labels, num_aln_files, elapsed_time, empty_alignment_files = find_and_process_alignments(
                        None, genbank_paths, args.cores, progress, None, None, args.analysis_type, args.all_alns_file
                    )
                else:  # gene homology
                    gene_homology_rows, num_aln_files, elapsed_time, empty_alignment_files, gene_count = find_and_process_alignments(
                        None, genbank_paths, args.cores, progress, None, None, args.analysis_type, args.all_alns_file
                    )
            else:
                # Existing directory-based processing
                if args.analysis_type == 'nucl':
                    simple_rows, detailed_rows, summary_df, np_matrix, df_matrix, contig_labels, num_aln_files, elapsed_time, empty_alignment_files = find_and_process_alignments(
                        args.alignments_dir, genbank_paths, args.cores, progress, expected_alignments, max_queue_size, args.analysis_type
                    )
                else:  # gene homology
                    gene_homology_rows, num_aln_files, elapsed_time, empty_alignment_files, gene_count = find_and_process_alignments(
                        args.alignments_dir, genbank_paths, args.cores, progress, expected_alignments, max_queue_size, args.analysis_type
                    )

            console.print(f"[green]Processed[/green] {num_aln_files} [green]alignment files in[/green] {elapsed_time:.2f} [green]seconds.[/green]")
            logging.info(f"[green]Processed[/green] {num_aln_files} [green]alignment files in[/green] {elapsed_time:.2f} [green]seconds.[/green]")

            if num_aln_files < expected_alignments:
                console.print(f"[bold red]Warning: Processed fewer alignment files than expected.[/bold red] [bold yellow]Expected:[/bold yellow] {expected_alignments}, [bold yellow]Processed:[/bold yellow] {num_aln_files}")
                logging.warning(f"[bold red]Warning: Processed fewer alignment files than expected.[/bold red] [bold yellow]Expected:[/bold yellow] {expected_alignments}, [bold yellow]Processed:[/bold yellow] {num_aln_files}")

            console.print("[yellow]Writing output files...[/yellow]")
            logging.info("Writing output files...")

            if args.analysis_type == 'nucl':
                write_output_files_nucl(
                    args.output_dir, simple_rows, detailed_rows, summary_df, np_matrix, df_matrix, contig_labels, progress, args.version, empty_alignment_files
                )
            else:  # gene homology
                write_output_files_gene(
                    args.output_dir, gene_homology_rows, progress, args.version, args.min_identity, args.identity_step, empty_alignment_files, gene_count
                )

            console.print(f"[green]Parsed results written to {args.output_dir} !![/green]")

            console.print(Panel.fit("[bold green]Homology Analysis Completed Successfully![/bold green]"))
            logging.info("[bold green]Homology Analysis Completed Successfully![/bold green]")

            final_log_path = os.path.join(args.output_dir, f'{args.logfile}.log')
            logging.info(f"Moving temp logfile to final logfile: {final_log_path}")
            shutil.move(temp_log_path, final_log_path)
            logging.info(f"Log file moved to: {final_log_path}")

    except Exception as e:
        print(f"An error occurred: {traceback.format_exc()}")

    finally:
        if os.path.exists(temp_log_path):
            os.unlink(temp_log_path)

if __name__ == "__main__":
    args = parse_arguments()
    if args.dale:
        args.cores = cpu_count() - 2
    elif args.cores:
        if args.cores <= 0:
            parser.error('ERROR! Please specify how many parallel cores to utilize via `--cpus` or use all cores with `--dale`')
            sys.exit(1)
        elif args.cores > cpu_count():
            parser.error("ERROR! Please don't specify more cores than are available. Run it in dale mode if you wanna redline this server!")
        elif args.cores == cpu_count():
            parser.error("ERROR! Please don't specify every single core. We need a few unburdened for OS stuff lest we risk memory corruption.\nRun it in dale mode if you wanna redline this server!")

    main(args)
