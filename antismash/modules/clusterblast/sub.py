# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Compares subsets of clusters to a reference data set of other subclusters
"""

import logging
import os
import resource
import time
from typing import Dict, List

from helperlibs.wrappers.io import TemporaryDirectory

from antismash.common import path
from antismash.common.secmet import Record
from antismash.config import ConfigType

from .core import (
    _SHIPPED_DATA_DIR,
    build_protein_cache,
    blastparse,
    get_core_gene_ids,
    load_clusterblast_database,
    run_blast,
    score_clusterblast_output,
    write_fastas_with_all_genes,
)
from .results import RegionResult, GeneralResults
from .data_structures import ProteinDB, ReferenceCluster


def _get_datafile_path(filename: str) -> str:
    """ A simple helper to get the full path of subclusterblast datafile """
    return os.path.join(_SHIPPED_DATA_DIR, 'sub', filename)


def prepare_sub_data(*, logging_only: bool = False) -> list[str]:
    """ Prepares the relevant data files.

        Arguments:
            logging_only: return error messages instead of raising exceptions

        Returns:
            a list of error messages
    """
    failure_messages = []
    try:
        build_protein_cache(_get_datafile_path("proteins.fasta"))
    except OSError as err:
        if logging_only:
            failure_messages.append(str(err))
        else:
            raise

    return failure_messages


def check_sub_prereqs(options: ConfigType) -> List[str]:
    """ Check if all required applications and datafiles are present.
        options is irrelevant here
    """
    _required_binaries = ['blastp', 'makeblastdb']

    _required_files = [
        'proteins.fasta',
        'proteins.fasta.phr',
        'proteins.fasta.pin',
        'proteins.fasta.psq',
        'clusters.txt'
    ]
    failure_messages = []
    for binary_name in _required_binaries:
        if binary_name not in options.executables:
            failure_messages.append(f"Failed to locate file: {binary_name!r}")

    for file_name in _required_files:
        if path.locate_file(_get_datafile_path(file_name)) is None:
            failure_messages.append(f"Failed to locate file: {file_name!r}")

    return failure_messages


def run_subclusterblast_on_record(record: Record, options: ConfigType) -> GeneralResults:
    """ Loads reference databases and performs subclusterblast analysis

        Arguments:
            record: the Record to analyse
            options: antismash Config

        Returns:
            a GeneralResults with results for each cluster in the record
    """
    logging.info('Running subcluster search')
    clusters, proteins = load_clusterblast_database("subclusterblast")
    return perform_subclusterblast(options, record, clusters, proteins)


def run_clusterblast_processes() -> None:
    """ Runs a single multi-threaded blastp over "input.fasta".

        blastp parallelises itself via ``-num_threads`` (see ``run_blast``), so
        no process-level fan-out is needed. This also works inside the daemonic
        Phase 2 worker processes, where a nested multiprocessing pool would
        silently fall back to running serially on a single core.

        Returns:
            None
    """
    database = _get_datafile_path('proteins.fasta')
    run_blast("input.fasta", database)


def read_clusterblast_output() -> str:
    """ Reads the blast output produced by run_clusterblast_processes.

        Returns:
            a string containing all blast run output
    """
    with open("input.out", encoding="utf-8") as handle:
        return handle.read()


def perform_subclusterblast(options: ConfigType, record: Record, clusters: Dict[str, ReferenceCluster],
                            proteins: ProteinDB) -> GeneralResults:
    """ Run BLAST on gene cluster proteins of each cluster, parse output and
        return result rankings for each cluster

        Arguments:
            options: antismash Config
            record: the Record to analyse
            clusters: a dictionary mapping reference cluster name to ReferenceCluster
            proteins: the protein database

        Returns:
            a GeneralResults instance storing results for all clusters in the
            record
    """
    results = GeneralResults(record.id, search_type="subclusterblast")
    blast_wall = 0.0
    blast_cpu = 0.0
    py_wall = 0.0
    py_cpu = 0.0
    with TemporaryDirectory(change=True):
        allcoregenes = get_core_gene_ids(record)
        for region in record.get_regions():
            # prepare and run blastp
            write_fastas_with_all_genes([region], "input.fasta")
            _w = time.time()
            _r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
            run_clusterblast_processes()
            _r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
            blast_wall += time.time() - _w
            blast_cpu += (_r1.ru_utime - _r0.ru_utime) + (_r1.ru_stime - _r0.ru_stime)
            # parse and score blastp results
            _w = time.time()
            _c0 = time.process_time()
            blastoutput = read_clusterblast_output()
            _, cluster_names_to_queries = blastparse(blastoutput, record,
                                                     min_seq_coverage=40,
                                                     min_perc_identity=45)
            ranking = score_clusterblast_output(clusters, allcoregenes, cluster_names_to_queries)
            # store results
            region_result = RegionResult(region, ranking, proteins, "subclusterblast")
            results.add_region_result(region_result, clusters, proteins)
            py_wall += time.time() - _w
            py_cpu += time.process_time() - _c0
    logging.info("CB-TIMING sub_blast: wall=%.3fs cpu=%.3fs cores=%.2f",
                 blast_wall, blast_cpu, blast_cpu / blast_wall if blast_wall > 0 else 0.0)
    logging.info("CB-TIMING sub_python: wall=%.3fs cpu=%.3fs cores=%.2f",
                 py_wall, py_cpu, py_cpu / py_wall if py_wall > 0 else 0.0)
    return results
