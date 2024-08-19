##############################################################
#                                                            #
# To whoever be reading this, peace be with ye and good luck #
#                                                            #
##############################################################
# Standard library imports
import argparse
import ast
import csv
import glob
import logging
import os
import re
import shutil
import sys
import tempfile
import time
import traceback
import pickle  # I need to add back the resume-from-crash feature I had in v2.
from collections import defaultdict, deque

# Stdlib multiprocessing imports
from concurrent.futures import (
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    as_completed,
)

from multiprocessing import cpu_count
from pathlib import Path # linter is angry but I've no clue why. Need to explicitly type the paths I guess?
from queue import Queue, Empty
from threading import Thread

# Third-party imports
import numpy as np
import pandas as pd
from Bio import SeqIO
from rich import box
from rich.console import Console
from rich.layout import Layout
from rich.live import Live
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.text import Text

# This is to get correct log printing behavior in our pane!
# instead of flushing the entire list once the capacity has been reached, it simply pops the oldest msg from the list.
class RollingBufferHandler(logging.Handler):
    """RollingBufferHandler

    Args:
        logging (logging.handler): give it a base logging.handler to wrap.
    This then creates a rolling buffer handler for our rich live progress panel to ensure that only the most recent
    messages are displayed!
    """
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
    """ConsolePanel
    Set up a custom __rich_console__ panel to use in our rich live progress panel.

    Args:
        Console (rich.console.Console): Give it the console handler to wrap. it then returns a corrected handler for the live
                                        view.
    """
    def __init__(self, *args, **kwargs):
        console_file = open(os.devnull, "w", encoding='utf-8')
        super().__init__(record=True, file=console_file, *args, **kwargs)

    def __rich_console__(
        self,
        console,
        options,
    ):  # stderr=False):
        texts = self.export_text(
            clear=False,
        ).split("\n")  # styles=True <- This fucks up the borders.
        for line in texts[-options.height :]:
            yield line


class LogPanel:
    def __init__(self, RollingBufferHandler):
        self.buffer_handler = RollingBufferHandler

    def __rich__(self) -> Panel:
        log_content = Text("\n").join(self.buffer_handler.get_logs())
        return Panel(log_content, box=box.ROUNDED, title="Recent Messages", expand=True)


# console = ConsolePanel()


def parse_genbank(path):
    """parse_genbank(path)

    Args:
        path (str or Path): the path to the genbank file of interest

    Returns:
        record: return SeqRecord for single contig in given file.
    {{TO-DO: Fix to generalize for multi-record files}}
    """
    try:
        records = list(SeqIO.parse(path, "genbank"))
        if not records:
            console.print(
                f"[bold yellow]Warning:[/bold yellow] Empty GenBank file: {path}"
            )
            return None
        return records[0]
    except Exception as e:
        console.print(
            f"[bold red]Error parsing GenBank file {path}:[/bold red] {str(e)}"
        )
        return None


def parse_all_genbanks(annotations_dir, progress):
    """parse_all_genbanks(annotations_dir, progres)

    Args:
        annotations_dir (str, Path): path to the directory containing the genbank files of interest
        progress (rich.progress.Progress): progress object to update global process bars in live view.

    Returns:
        dict:
            key: genbank_id,
            value: list of paths for each given genbank
    """
    genbank_files = glob.glob(f"{annotations_dir}/*.gbff")
    genbank_paths = {}
    total_files = len(genbank_files)
    task = progress.add_task(
        "[cyan]Parsing genbanks...", total=total_files, expand=True
    )
    for i, file in enumerate(genbank_files):
        assembly_contig_id = os.path.basename(file).split(".")[0]
        genbank_paths[assembly_contig_id] = file
        progress.update(
            task,
            advance=1,
            description=f"[cyan]Indexing GenBank files... ({i+1}/{total_files})",
            expand=True,
        )
    return genbank_paths


def parse_alignment(path):
    """parse_alignment(path)

    Args:
        path (str, Path): path to the alignment file of interest

    Returns:
        dict:
            key: alignment file columns,
            values: corresponding values
    """
    with open(path, "r", encoding='utf-8') as infile:
        lines = infile.readlines()
    keys = lines[0].strip().split("\t")
    return [dict(zip(keys, line.strip().split("\t"))) for line in lines[1:]]

def is_within_range(feature, start, end):
    """Check to see if a genomic feature is within the interval we are looking at"""
    return int(feature.location.start) >= start and int(feature.location.end) <= end

def get_features_from_range(record, start, end):
    """Look at all features that are present within start:end for our given record."""
    return [
        feature
        for feature in record.features
        if is_within_range(feature, start, end) and feature.type == "CDS"
    ]

def simplify_genes_for_contig(features):
    """simplify the genes returned for each contig given a list of features.
    This can probably be deleted at some point"""
    return [
        (
            feature.qualifiers["locus_tag"][0].strip("'[]"),
            " ".join(feature.qualifiers["product"]),
        )
        for feature in features
    ]


def get_coverage(length, coords):
    """get_coverage(length, coords)
    Given the length of a contig and the coordinates of alignments, determine percent coverage, identities, and length of homology.

    Args:
        length (int): length of the source contig
        coords (list): list of tuples containing (start, end) for each region of alignment.

    Returns:
        percent_coverage(float): percentage of source contig covered by alignments.
        identities (int): number of identities between source and target contigs.
        total_homology(int): integer value of total length of shared homology.
    """
    covered_positions = np.zeros(length, dtype=bool)
    identities = 0
    for start, end, ident in coords:
        if start > end:
            start, end = end, start
        covered_positions[start : end + 1] = True
        identities += (end - start + 1) * float(1 / ident)
    total_homology = np.sum(covered_positions)
    percent_coverage = total_homology / length * 100
    return percent_coverage, identities, total_homology


class CursedDict(defaultdict):
    """
    may the divine forgive me for such an egregious hack.
    Basically to keep from refactoring everything via a custom factory we just
    commit this atrocity and instantiate two hard coded keys and then treat the
    rest as a defaultdict(list).
    This is going to break many thing probably but I hope it fixes at least two problems.
    """

    def __init__(self):
        super().__init__(list)  # super().__scuffed__
        self["assembly_id"] = ""
        self["contig_length"] = int()

    #def __getitem__(self, key):  # no hardcoding required? this is unnecessary per the all knowing linter.
        # if key in ("assembly_id", "contig_length"):
        #    return super().__getitem__(key)
        #return super().__getitem__(key)

    def __setitem__(self, key, value):
        # this is so bad lmao
        if key == "assembly_id" and not isinstance(value, str):
            raise ValueError("assembly_id must be a string")
        elif key == "contig_length" and not isinstance(value, int):
            raise ValueError("contig_length must be an integer")
        super().__setitem__(key, value)


def get_homologies(alignments, asm1_record, asm2_record):
    """get_homologies

    Args:
        alignments (dict): dictionary generated by parse_alignments
        asm1_record (SeqRecord): SeqRecord of parsed genbank for assembly1
        asm2_record (SeqRecord): SeqRecord of parsed genbank for assembly2

    Returns:
        genes_dict ( defaultdict(CursedDict) ): dictionary containing alignments between source(asm1) and target(asm2) where each alignment is a dictionary
        detailing the specific alignment at a set of coordinates, top level key is asm1, subkeys for asm1 are 'contig_length' and 'assembly_id'
            # DO NOT CHANGE THIS IT WILL BREAK DOWNSTREAM PARSING.
    """
    genes_dict = defaultdict(CursedDict)
    for alignment in alignments:
        asm1_id = alignment["QUERY_ID"]
        asm1_name = alignment["QUERY_NAME"]
        asm1_start = int(alignment["QUERY_START"])
        asm1_end = int(alignment["QUERY_END"])
        asm1_length = int(len(asm1_record.seq))
        asm2_id = alignment["REF_ID"]
        asm2_name = alignment["REF_NAME"]
        asm2_start = int(alignment["REF_START"])
        asm2_end = int(alignment["REF_END"])
        asm2_length = int(len(asm2_record.seq))
        asm2_identity = float(alignment["IDENTITY"])

        asm1_features = get_features_from_range(asm1_record, asm1_start, asm1_end)
        asm2_features = get_features_from_range(asm2_record, asm2_start, asm2_end)
        asm1_genes = simplify_genes_for_contig(asm1_features)
        asm2_genes = simplify_genes_for_contig(asm2_features)

        genes_dict[asm1_id]["contig_length"] = asm1_length
        genes_dict[asm1_id]["assembly_id"] = alignment["QUERY_ID"]

        if asm2_id not in genes_dict[asm1_id]:
            genes_dict[asm1_id][asm2_id] = []

        genes_dict[asm1_id][asm2_id].append(
            {
                "asm_aln": {
                    "asm_name": asm1_name,
                    "start": asm1_start,
                    "end": asm1_end,
                    "aln_length": int(alignment["QUERY_LENGTH"]),
                    "features": asm1_genes,
                },
                "ref_aln": {
                    "ref_id": asm2_id,
                    "ref_name": asm2_name,
                    "start": asm2_start,
                    "end": asm2_end,
                    "aln_length": int(alignment["REF_LENGTH"]),
                    "aln_identity": asm2_identity,
                    "features": asm2_genes,
                    "ref_length": asm2_length,
                },
            }
        )
    return genes_dict


def make_table_for_asm(genes_dict):
    """using the dictionary of homologies generated by get_homologies(), output lines for the (detailed_rows) list"""
    lines = []
    for asm1_id, asm1_data in genes_dict.items():
        assembly1_id = asm1_data["assembly_id"]
        asm1_len = asm1_data["contig_length"]
        for asm2_id, alignments in asm1_data.items():
            if asm2_id not in ["contig_length", "assembly_id"]:
                asm1_aln_coords = [
                    (
                        cov["asm_aln"]["start"],
                        cov["asm_aln"]["end"],
                        cov["ref_aln"]["aln_identity"],
                    )
                    for cov in alignments
                ]
                asm2_aln_coords = [
                    (
                        cov["ref_aln"]["start"],
                        cov["ref_aln"]["end"],
                        cov["ref_aln"]["aln_identity"],
                    )
                    for cov in alignments
                ]
                asm2_len = alignments[0]["ref_aln"]["ref_length"]
                assembly2_id = alignments[0]["ref_aln"]["ref_id"]
                asm2_genes = [
                    gene for item in alignments for gene in item["ref_aln"]["features"]
                ]
                asm1_genes = [
                    gene for item in alignments for gene in item["asm_aln"]["features"]
                ]
                asm1_cov, asm1_identities, asm1_total_homology = get_coverage(
                    asm1_len, asm1_aln_coords
                )
                asm2_cov, asm2_identities, asm2_total_homology = get_coverage(
                    asm2_len, asm2_aln_coords
                )
                lines.append(
                    [
                        assembly1_id,
                        asm1_id,
                        asm1_len,
                        assembly2_id,
                        asm2_id,
                        asm2_len,
                        asm1_cov,
                        asm1_total_homology,
                        asm1_identities,
                        asm2_cov,
                        asm2_total_homology,
                        asm2_identities,
                        len(asm1_genes),
                        len(asm2_genes),
                    ]
                )
    return lines


def get_rows(alignments, genes_dict):
    """using our dictionary of alignments parsed from our alignments file via get_alignments"""
    rows = []
    for contig, contig_data in genes_dict.items():
        contig_len = contig_data["contig_length"]
        assembly_id = contig_data["assembly_id"]
        total_contig_cov = 0
        contig_alignments = {}

        # Handle the case where contig_len appears as a list
        if isinstance(contig_len, list):
            contig_len = contig_len[
                0
            ]  # Take the first element if it exists, otherwise use 0
        contig_total_homology = 0
        for aln, alignments in contig_data.items():
            if aln not in ["contig_length", "assembly_id"]:
                asm_coords = [
                    (
                        cov["asm_aln"]["start"],
                        cov["asm_aln"]["end"],
                        cov["ref_aln"]["aln_identity"],
                    )
                    for cov in alignments
                ]
                ref_coords = [
                    (
                        cov["ref_aln"]["start"],
                        cov["ref_aln"]["end"],
                        cov["ref_aln"]["aln_identity"],
                    )
                    for cov in alignments
                ]
                ref_len = alignments[0]["ref_aln"]["ref_length"]
                _, _, contig_total_homology = get_coverage(contig_len, asm_coords)
                #_, _, ref_total_homology = get_coverage(ref_len, ref_coords)
                contig_alignments[aln] += int(contig_total_homology)
                total_contig_cov += contig_total_homology

        # Only add the row if there are alignments and the contig has non-zero length
        if contig_alignments and contig_len > 0:
            rows.append(
                [
                    assembly_id,
                    contig,
                    int(contig_len),  # Ensure this is an integer
                    int(total_contig_cov),  # Ensure this is an integer
                    [
                        f"({k}:{v})" for k, v in contig_alignments.items()
                    ],  # Format alignments as a list of strings
                ]
            )
    return rows


def check_tsv_content(file_path):
    """Make sure a tsv file isn't empty and has more than just a header row."""
    with open(file_path, "r", encoding='utf-8') as file:
        lines = file.readlines()
        return True if len(lines) > 1 else False


def process_alignment(alignment_file, genbank_paths):
    """given an alignment file, process the alignments within and use the provided genbanks to also pull required information
    about the two genomes being compared."""
    asm1_id, asm2_id = os.path.basename(os.path.dirname(alignment_file)).split("_vs_")

    asm1_record = parse_genbank(genbank_paths[asm1_id])
    asm2_record = parse_genbank(genbank_paths[asm2_id])

    if not check_tsv_content(alignment_file):
        # File is empty, but we still want to return a structured result
        simple_row = [asm1_id, asm1_id, len(asm1_record.seq), "0", []]
        detailed_row = [
            asm1_id,
            asm1_id,
            len(asm1_record.seq),
            asm2_id,
            asm2_id,
            len(asm2_record.seq),
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
        return (
            True,
            f"NO HOMOLOGY BETWEEN {asm1_id} and {asm2_id}!",
            [simple_row],
            [detailed_row],
        )

    alignments = parse_alignment(alignment_file)
    homologies = get_homologies(alignments, asm1_record, asm2_record)
    detailed_rows = make_table_for_asm(homologies)
    simple_rows = get_rows(alignments, homologies)

    return (
        True,
        f"Parsed homology between {asm1_id} and {asm2_id}!",
        simple_rows,
        detailed_rows,
    )

def discover_files(directory, queue):
    """given a directory and a queue (for multiprocessing), walk through it and find all paths for the file specified below
    returns nothing as all discovered paths are put into the queue used by other processes."""
    discovered_count = 0
    logging.debug("discover_files: Beginning walk through %s", directory)
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename == "align_coords.tsv":
                full_path = os.path.join(dirpath, filename)
                queue.put(full_path)
                discovered_count += 1
    # Signal that discovery is complete
    logging.debug("discover_files: Walk finished. Inserting None into queue.")
    queue.put(None)
    logging.debug("discover_files: Discovery complete. Found %s files.", discovered_count)

def worker(file_list, genbank_paths):
    """Worker function to execute alignments given a list of files and list of genbank paths.
    This is executed in parallel as a big multiprocessing pool party.
    This worker spins up as many threads as the system will allow it to create. Each thread is given a file to process.
    all threads then return results and the total list of results is returned to the multiprocessing executor.
    These results are then accessed through as_completed(futures) on the main executor process."""
    results = []
    with ThreadPoolExecutor() as thread_executor:
        futures = [
            thread_executor.submit(process_alignment, file_path, genbank_paths)
            for file_path in file_list
        ]
        for future in as_completed(futures):
            results.append(future.result())
    return results

def find_and_process_alignments(
    directory, genbank_paths, num_cpus, progress, expected_count, max_queue_size
):
    """this performs multiprocessing and multithreading to locate paths to alignment files, process them in parallel, and then
    output results in lists, matrices, and dataframes.

    a separate thread for discover_files is spun up prior to beginning multiprocessing execution.
    Once a specified queue size has been reached, jobs are then scheduled and distributed to each parallel worker,
    on each worker, the list of files and genbank paths is further divided into their own individual threads.

    once all files have been discovered and added to the queue, the discovery thread is merged back with the main thread
    and the remaining processing is run to completion.

    ( I am proud of this but it could do with much more optimization, I'm positive of this.)

    Max queue size was -un-implemented as process count * 10000 just works the best
    but this could be optimized and reimplemented in the future.
    """
    start_time = time.time()
    simple_rows = [
        [
            "assembly_id",
            "contig_id",
            "contig_len",
            "total_contig_homology",
            "list_of_alignments(asm2:homology_length)",
        ]
    ]
    detailed_rows = [
        [
            "assembly1_id",
            "asm1_id",
            "asm1_len",
            "assembly2_id",
            "asm2_id",
            "asm2_len",
            "asm1_cov",
            "asm1_total_homology",
            "asm1_identities",
            "asm2_cov",
            "asm2_total_homology",
            "asm2_identities",
            "asm1_genes",
            "asm2_genes",
        ]
    ]

    logging.debug("Beginning parallel execution")

    logging.debug("Beginning file discovery, initializing file Queue() and discovery Thread()")
    discovery_task = progress.add_task(
        "[cyan]Discovering alignment files...", total=None, expand=True
    )
    processing_task = progress.add_task(
        "[cyan]Processing alignment files...", total=expected_count, expand=True
    )

    file_queue = Queue()
    discovery_thread = Thread(target=discover_files, args=(directory, file_queue)) # start discovery thread, feed it the directory and file queue.

    logging.debug("Starting discovery thread!")
    discovery_thread.start()

    discovered_count = 0
    processed_count = 0

    with ProcessPoolExecutor(max_workers=num_cpus) as process_executor: # start the pool party.
        futures = []
        discovery_complete = False
        all_files_processed = False
        file_buffer = []

        while not discovery_complete or futures or file_buffer:
            while not discovery_complete and len(file_buffer) < num_cpus * 10000: # see docstring for future optimization.
                try:
                    aln_path = file_queue.get(timeout=1)
                    if aln_path is None:
                        # branch 1
                        discovery_complete = True
                        logging.debug("`None` path reached! Discovery is finished. Waiting for discovery thread to complete. Merging thread with main!")
                        discovery_thread.join()
                        logging.debug("Discovery thread successfully merged!")
                        progress.update(
                            discovery_task,
                            completed=discovered_count,
                            total=discovered_count,
                            description=f"[bold green]File discovery complete! Found[/bold green] {discovered_count} [bold green]alignments.[/bold green]",
                        )
                        break

                    file_buffer.append(aln_path)
                    discovered_count += 1
                except Empty:
                    # empty branch
                    if not discovery_thread.is_alive():
                        discovery_complete = True
                        logging.debug(
                            "discovery thread is dead and queue is empty. merging the corpse with main thread."
                        )
                        discovery_thread.join()
                        logging.debug("discovery thread merged with main.")
                        progress.update(
                            discovery_task,
                            completed=discovered_count,
                            total=discovered_count,
                            description=f"[bold green]File discovery complete! Found[/bold green] {discovered_count} [bold green]alignments.[/bold green]",
                        )
                        break

            while file_buffer and len(futures) < num_cpus:
                chunk_size = len(file_buffer) // num_cpus + 1
                file_chunk = file_buffer[:chunk_size]
                file_buffer = file_buffer[chunk_size:]
                futures.append(
                    process_executor.submit(worker, file_chunk, genbank_paths)
                )

            for future in as_completed(futures):
                futures.remove(future)
                try:
                    results = future.result()
                    for result in results:
                        simple_rows.extend(result[2])
                        detailed_rows.extend(result[3])

                    processed_count += len(results)
                    progress.update(
                        processing_task,
                        advance=len(results),
                        description=f"[cyan]Processing alignment files... ({processed_count} processed)",
                    )
                except Exception:
                    logging.error(
                        "Error during file processing:\n %s", traceback.format_exc()
                    )

                if processed_count >= expected_count and not all_files_processed:
                    logging.info("Reached expected count of %s, setting all_files_processed flag to True!", expected_count)
                    all_files_processed = True
                    progress.update(processing_task, total=discovered_count)

    elapsed_time = time.time() - start_time
    logging.info(
        "Parallel processing is complete! Total Runtime: %s seconds. Files discovered: %i, Files processed: %i",
        f'{elapsed_time:.2f}', discovered_count, processed_count
    )

    logging.info("Total rows in simple_rows: %s", (len(simple_rows) - 1))  # Subtract 1 for the header row

    # Aggregate alignments
    logging.info("creating dataframe from simple_rows")
    df = pd.DataFrame(simple_rows[1:], columns=simple_rows[0])

    logging.info("Total rows in DataFrame: %s", len(df))

    logging.info("aggregating alignments by contig")
    grouped_data = aggregate_alignments(df)

    logging.info("Total contigs in grouped_data: %s", len(grouped_data))

    # Create summary DataFrame
    logging.info("creating summary dataframe of aggregated alignments")
    summary_df = create_summary_dataframe(grouped_data)

    logging.info("Total rows in summary DataFrame: %s", len(summary_df))

    # Create homology matrix
    logging.info("creating matrices from aggregated alignments!")
    np_matrix, df_matrix, contig_labels = create_homology_matrix(grouped_data)

    # Ensure unscuffed matrix diagonal of the np_matrix, df_matrix, and throw tantrum if broken!
    zero_contigs_np = check_matrix_diagonal(np_matrix, matrix_type='np', contig_labels=contig_labels)
    zero_contigs_df = check_matrix_diagonal(df_matrix, matrix_type='df')

    if zero_contigs_np:
        logging.error("The following np_contigs have zero length (diagonal value): %s", ', '.join(zero_contigs_np))
    if zero_contigs_df:
        logging.error("The following df_contigs have zero length (diagonal value): %s", ', '.join(zero_contigs_df))
    else:
        logging.info("All contigs have non-zero lengths (diagonal values).")

    return (simple_rows, detailed_rows, summary_df, np_matrix, df_matrix, contig_labels, processed_count, elapsed_time)

def aggregate_alignments(df):
    """using our dataframe created from simple_rows, aggregate all of the alignments for each contig
    and return a dictionary of grouped alignments by each contig."""
    grouped = {}
    all_contigs = set()

    for _, row in df.iterrows():
        source_contig = row["contig_id"]
        all_contigs.add(source_contig)

        if source_contig not in grouped:
            grouped[source_contig] = {
                "source_len": row["contig_len"],
                "total_contig_homology": row["total_contig_homology"],
                "alignments": {},
            }

        alignments = row["list_of_alignments(asm2:homology_length)"]
        if isinstance(alignments, str):
            alignments = ast.literal_eval(alignments)  # eval() is dangerous, use ast.literal_eval() instead for ""Safety""

        for alignment in alignments:
            match = re.match(r"\(([^:]+):(\d+)\)", alignment)
            if match:
                target_contig, homology_len = match.groups()
                homology_len = int(homology_len)
                all_contigs.add(target_contig)
                grouped[source_contig]["alignments"][target_contig] = homology_len

                # Add target contig to grouped if it's not there
                if target_contig not in grouped:
                    grouped[target_contig] = {
                        "source_len": None,  # We don't know the length yet
                        "total_contig_homology": None,
                        "alignments": {},
                    }

    # Ensure all contigs are in the grouped dictionary
    for contig in all_contigs:
        if contig not in grouped:
            grouped[contig] = {
                "source_len": None,
                "total_contig_homology": None,
                "alignments": {},
            }

    return grouped

def create_summary_dataframe(grouped_data):
    """using our grouped data dictionary, turn that into a summary dataframe"""
    summary_data = []
    logging.debug("beginning summary df construction!")
    for source_contig, data in grouped_data.items():
        summary_data.append(
            {
                "source_contig": source_contig,
                "source_len": data["source_len"],
                "total_contig_homology": data["total_contig_homology"],
                "num_alignments": len(data["alignments"]),
                "target_contigs": ", ".join(data["alignments"].keys()),
                "homology_lengths": ", ".join(map(str, data["alignments"].values())),
            }
        )
    logging.debug("finished summary df construction!")
    return pd.DataFrame(summary_data)

def create_homology_matrix(detailed_rows):
    """Creates a homology matrix from the detailed rows.
    this function has removed at least a year off of my total lifespan"""
    # Extract all unique contigs
    all_contigs = set()
    contig_lengths = {}
    for row in detailed_rows[1:]:  # Skip header row
        asm1_id, asm1_len, asm2_id, asm2_len = row[1], int(row[2]), row[4], int(row[5])
        all_contigs.add(asm1_id)
        all_contigs.add(asm2_id)
        contig_lengths[asm1_id] = asm1_len
        contig_lengths[asm2_id] = asm2_len

    all_contigs = sorted(list(all_contigs))
    n = len(all_contigs)

    # Create a mapping of contig names to matrix indices
    contig_to_index = {contig: i for i, contig in enumerate(all_contigs)}

    # Initialize the matrix with zeros
    np_matrix = np.zeros((n, n), dtype=np.int64)

    # Fill the diagonal with contig lengths
    for i, contig in enumerate(all_contigs):
        np_matrix[i, i] = contig_lengths[contig]

    # Fill the matrix with homology data
    for row in detailed_rows[1:]:  # Skip header row
        asm1_id, asm2_id, asm1_total_homology = row[1], row[4], int(row[7])
        i, j = contig_to_index[asm1_id], contig_to_index[asm2_id]
        np_matrix[i, j] = asm1_total_homology
        np_matrix[j, i] = asm1_total_homology  # Ensure symmetry

    # Create pandas DataFrame
    df_matrix = pd.DataFrame(np_matrix, index=all_contigs, columns=all_contigs)

    return np_matrix, df_matrix, all_contigs

def check_matrix_diagonal(matrix, matrix_type, **kwargs):
    """
    Check the diagonal cells of the given homology matrix to ensure none are zero.
    Can be used for the numpy matrix but MUST give labels as kwarg.
    Can be used for the pandas matrix as well.

    Args:
    np_matrix (numpy.ndarray or pandas.DataFrame): The matrix representing the homology data.
    matrix_type (str): either 'np' or 'df' for numpy array or pandas dataframe.

    Kwargs:
    contig_labels (list): List of contig labels corresponding to matrix rows/columns. Required if matrix_type='np'.

    Returns:
    list: List of contig labels with zero diagonal values, if any.
    """
    if matrix_type == 'np':
        if 'contig_labels' not in kwargs:
            raise ValueError("contig_labels must be provided as a kwarg when matrix_type is 'np'")

        contig_labels = kwargs['contig_labels']

        if len(contig_labels) != matrix.shape[0]:
            raise ValueError("Length of contig_labels must match the matrix dimensions")

        zero_diagonals = [contig for i, contig in enumerate(contig_labels) if matrix[i, i] == 0]
        return zero_diagonals

    elif matrix_type in ['pd', 'df']:
        zero_diagonals = [contig for contig in matrix.index if matrix.loc[contig, contig] == 0]
        return zero_diagonals

    else:
        raise ValueError("matrix_type must be either 'np' or 'df'")

def write_output_files(
    output_dir,
    simple_rows,
    detailed_rows,
    summary_df,
    np_matrix,
    df_matrix,
    contig_labels,
    progress,
    version,
):
    """self explanatory function. give it the output directory and the various args as calculated then output within that given dir."""
    logging.info("Writing rows to files. Please stand by.")
    task = progress.add_task("[cyan]Writing output files...", total=6)

    # Write simple_rows
    logging.info("Writing simple rows!")
    with open(f"{output_dir}/ava_homology_simple_{version}.tsv", "w", newline="", encoding='utf-8') as f:
        writer = csv.writer(f, delimiter="\t")
        for row in simple_rows:
            writer.writerow(row)
    progress.update(task, advance=1)
    logging.info("Simple_rows file has been written!")

    # Write detailed_rows
    logging.info("Writing detailed rows!")
    with open(f"{output_dir}/ava_homology_detailed_{version}.tsv", "w", newline="", encoding='utf-8',) as f:
        writer = csv.writer(f, delimiter="\t")
        for row in detailed_rows:
            writer.writerow(row)
    progress.update(task, advance=1)
    logging.info("Detailed_rows file has been written!")

    # Write summary DataFrame
    logging.info("Writing summary dataframe!")
    summary_df.to_csv(
        f"{output_dir}/ava_homology_summary_{version}.tsv", sep="\t", index=False, encoding='utf-8',
    )
    progress.update(task, advance=1)
    logging.info("Summary DataFrame has been written!")

    # Write numpy matrix
    logging.debug("Writing numpy matrix!")
    np.savetxt(
        f"{output_dir}/ava_homology_matrix_{version}.csv", np_matrix, delimiter=","
    )
    progress.update(task, advance=1)
    logging.info("Numpy matrix has been written!")

    # Write pandas DataFrame matrix
    logging.debug("Writing pandas matrix!")
    df_matrix.to_csv(f"{output_dir}/ava_homology_matrix_labeled_{version}.csv")
    progress.update(task, advance=1)
    logging.info("Pandas matrix has been written!")

    # Write contig labels
    logging.debug("writing contig labels (for numpy row/cols)")
    with open(f"{output_dir}/ava_homology_np_matrix_labels_{version}.txt", "w", encoding='utf-8',) as f:
        for label in contig_labels:
            f.write(f"{label}\n")
    progress.update(task, advance=1)
    logging.info("Contig labels have been written!")

    logging.info("All output files have been written!")

def setup_logging(log_file):
    """Set up logging to both file and console."""
    logger = logging.getLogger()
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )
    return logger

def create_layout(console, progress, args_str, buffering_handler):
    """Creates a rich cli layout with quad panes as follows:
    top_right: console output (stdout & stderr)
    top_left: logging output (using a rolling buffer handler) ((still logs to file))
    bottom_right: Progress bars for each task
    bottom_left: Arguments provided for execution.
    """
    log_panel = LogPanel(buffering_handler)

    layout = Layout()
    layout.split(Layout(name="top", ratio=3), Layout(name="bottom", ratio=1))
    layout["top"].split_row(
        Layout(
            Panel(console, title="Console Output", box=box.HEAVY, expand=True),
            name="top_left",
            ratio=1,
        ),
        Layout(
            Panel(log_panel, title="DEBUG", box=box.HEAVY, expand=True),
            name="top_right",
            ratio=2,
        ),
    )
    layout["bottom"].split_row(
        Layout(
            Panel(
                progress,
                title="Task Progress",
                box=box.HEAVY,
            ),
            name="bottom_right",
            ratio=2,
        ),
        Layout(
            Panel.fit(args_str, title="Arguments", box=box.HEAVY),
            name="bottom_left",
            ratio=1,
        ),
    )
    return layout


def parse_arguments():
    """Parse argparse cli options."""
    parser = argparse.ArgumentParser(description="Homology Analysis")
    parser.add_argument(
        "--annotations_dir", required=True, help="Directory containing GenBank files"
    )
    parser.add_argument(
        "--alignments_dir", required=True, help="Directory containing alignment files"
    )
    parser.add_argument(
        "--output_dir", required=True, help="Directory to store output files"
    )
    parser.add_argument("--logfile", required=True, help="Path to log file")
    parser.add_argument(
        "--cores", type=int, help="Number of CPU cores to use", required=False
    )
    parser.add_argument(
        "--dale",
        action="store_true",
        help="Redline this server by using all cores (less two for headroom)",
        required=False,
    )
    parser.add_argument(
        "--version", type=str, help="Which version is this?", required=True
    )
    return parser.parse_args()


def main(args):
    """run entire analysis using the provided cli args parsed by parse_arguments()"""
    console = ConsolePanel()

    # Create a temporary file for logging
    temp_log_file = tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".log")
    temp_log_path = temp_log_file.name
    temp_log_file.close()  # Close the file but don't delete it
    try:
        main_logger = setup_logging(temp_log_path)  # set up to temp logfile.
        main_logger.setLevel(logging.DEBUG)
        logging.debug("initialized logger: main_logger")

        # Instantiate the buffering handler
        LOG_BUFFER_MAX_MSGS = 16
        logging.debug("max buffer length: %i", LOG_BUFFER_MAX_MSGS)
        logging.debug("setting up log formatter for pane view")
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
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
        args_dict = vars(args)  # Convert args to dict
        args_str = "\n".join(
            f"[bold green]{k}[/bold green]: {v}" for k, v in args_dict.items()
        )
        logging.debug(
            "creating path for subdirectory:%s in output:%s", args.version, args.output_dir
        )
        args.output_dir = os.path.join(
            args.output_dir, args.version
        )  # make a version subdirectory for some order in this chaos
        logging.debug("new output path: %s", args.output_dir)

        layout = create_layout(console, progress, args_str, buffering_handler)

        logging.debug("Starting live view")
        with Live(layout, refresh_per_second=10) as live:
            logging.info("Using the following arguments for this run\n%s\n", args_str)

            logging.info("Creating output directory: %s", args.output_dir)
            console.print(f"Creating output directory: {args.output_dir}")
            os.makedirs(args.output_dir, exist_ok=True)
            logging.info("Output directory successfully created")
            console.print("[green]Output directory successfully created[/green]")

            console.print(
                Panel.fit(
                    "[bold cyan]Starting Homology Analysis[/bold cyan]",
                    highlight=True,
                    box=box.ROUNDED,
                )
            )
            logging.info("Starting Homology Analysis")

            logging.info("Indexing all GenBank files...")
            genbank_paths = parse_all_genbanks(args.annotations_dir, progress)

            if not genbank_paths:
                console.print(
                    "[bold red]Error: No valid GenBank files found. Exiting.[/bold red]"
                )
                logging.error(
                    "Error: No valid GenBank files found. Exiting."
                )
                sys.exit(1)

            n_genbanks = len(genbank_paths)
            console.print(
                f"[green]Successfully indexed[/green] [green]{n_genbanks}[/green] [green]GenBank files.[/green]"
            )
            logging.info("Successfully indexed %s GenBank files.", n_genbanks)

            expected_alignments = (n_genbanks * (n_genbanks - 1)) // 2 # N*(N-1)/2 alignments expected
            console.print(
                f"[yellow]Expected number of alignment files:[/yellow] {expected_alignments}"
            )
            logging.info("Expected number of alignment files: %s", expected_alignments)

            console.print(
                "[bold blue]Finding and processing alignment files...[/bold blue]"
            )
            logging.info("Finding and processing alignment files...")
            max_queue_size = 10000  # cli arg? need to optimize for core count and maybe functionalize that.
            (
                simple_rows,
                detailed_rows,
                summary_df,
                np_matrix,
                df_matrix,
                contig_labels,
                num_aln_files,
                elapsed_time,
            ) = find_and_process_alignments(
                args.alignments_dir,
                genbank_paths,
                args.cores,
                progress,
                expected_alignments,
                max_queue_size,
            )

            console.print(f"[green]Processed[/green] {num_aln_files} [green]alignment files in[/green] {elapsed_time:.2f} [green]seconds.[/green]")
            logging.info("Processed %i alignment files in %s seconds.", num_aln_files, f'{elapsed_time:.2f}')

            if num_aln_files < expected_alignments:
                console.print(f"[bold red]Warning: Processed fewer alignment files than expected.[/bold red] [bold yellow]Expected:[/bold yellow] {expected_alignments}, [bold yellow]Processed:[/bold yellow] {num_aln_files}")
                logging.warning(
                    "Warning: Processed fewer alignment files than expected.Expected: %i, Processed: %i",
                    expected_alignments, num_aln_files
                )

            console.print("[yellow]Writing output files...[/yellow]")
            logging.info("Writing output files...")

            write_output_files(
                args.output_dir,
                simple_rows,
                detailed_rows,
                summary_df,
                np_matrix,
                df_matrix,
                contig_labels,
                progress,
                args.version,
            )

            console.print(
                f"[green]Parsed results written to {args.output_dir} !![/green]"
            )

            console.print(
                Panel.fit(
                    "[bold green]Homology Analysis Completed Successfully![/bold green]"
                )
            )
            logging.info(
                "Homology Analysis Completed Successfully!"
            )

            # Move the temporary log file to the output directory
            final_log_path = os.path.join(args.output_dir, f"{args.logfile}.log")
            logging.info("Moving temp logfile to final logfile:", final_log_path)
            shutil.move(temp_log_path, final_log_path)
            logging.info("Log file moved to: ", final_log_path)

    except Exception:
        print(f"An error occurred: {traceback.format_exc()}")

    finally:
        # Ensure the temporary file is removed if it still exists
        if os.path.exists(temp_log_path):
            os.unlink(temp_log_path)


if __name__ == "__main__":
    args = parse_arguments()
    if args.dale: # do it for Dale
        args.cores = cpu_count() - 2  # give all cores less two for some headroom :)
    elif args.cores:
        if args.cores <= 0:
            parser.error(
                "ERROR! Please specify how many parallel cores to utilize via `--cpus` or use all cores with `--dale`"
            )
            sys.exit(1)
        elif args.cores > cpu_count():
            parser.error(
                "ERROR! Please don't specify more cores than are available. Run it in dale mode if you wanna redline this server!"
            )
        elif args.cores == cpu_count():
            parser.error(
                "ERROR! Please don't specify every single core. We need a few unburdened for OS stuff lest we risk memory corruption.\nRun it in dale mode if you wanna redline this server!"
            )

    main(args)
