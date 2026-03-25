# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running hmmsearch.
"""

import logging
import multiprocessing
import os
import re
from typing import List

from helperlibs.wrappers.io import TemporaryDirectory

from .base import execute, get_config, get_effective_cpus, SearchIO, parallel_function


def run_hmmsearch(query_hmmfile: str, target_sequence: str, use_tempfile: bool = False,
                  force_single_core: bool = False
                  ) -> List[SearchIO._model.query.QueryResult]:  # pylint: disable=protected-access
    """ Run hmmsearch on a HMM file and a fasta input

        Arguments:
            query_hmmfile: the path to the HMM file
            target_sequence: the fasta input to search as a string
            use_tempfile: if True, a tempfile will be written for the fasta input
                          instead of piping
            force_single_core: if True, only use 1 CPU core regardless of config

        Returns:
            a list of hmmsearch results as parsed by SearchIO
    """
    config = get_config()

    # Check if multithreading is disabled in the config
    multithreading_disabled = (config.get('hmmer3') and 'multithreading' in config.hmmer3
                               and not config.hmmer3.multithreading)

    effective_cpus = get_effective_cpus()

    # If multithreading is disabled, only 1 CPU is available, or force_single_core is True,
    # use the non-parallel version
    if multithreading_disabled or effective_cpus == 1 or force_single_core or multiprocessing.current_process().daemon:
        cpu_count = 1 if force_single_core else effective_cpus

        command = [config.executables.hmmsearch, "--cpu", str(cpu_count),
                   "-o", os.devnull,  # throw away the verbose output
                   "--domtblout", "result.domtab",
                   query_hmmfile]

        # Remove CPU flag if multithreading is disabled
        if multithreading_disabled:
            command = command[0:1] + command[3:]

        with TemporaryDirectory(change=True):
            try:
                if use_tempfile:
                    with open("input.fa", "w", encoding="utf-8") as handle:
                        handle.write(target_sequence)
                    command.append("input.fa")
                    run_result = execute(command)
                else:
                    command.append('-')
                    run_result = execute(command, stdin=target_sequence)
            except OSError:
                return []
            if not run_result.successful():
                logging.error('hmmsearch returned %d: %s while searching %s',
                              run_result.return_code, run_result.stderr, query_hmmfile)
                raise RuntimeError("Running hmmsearch failed.")
            return list(SearchIO.parse("result.domtab", 'hmmsearch3-domtab'))
    else:
        # Use the parallel implementation
        return run_hmmsearch_parallel(query_hmmfile, target_sequence, use_tempfile)


def split_fasta(fasta_text: str) -> List[str]:
    """ Split a multi-record FASTA string into individual records

        Arguments:
            fasta_text: a string containing multiple FASTA records

        Returns:
            a list of strings, each containing a single FASTA record
    """
    if not fasta_text:
        return []

    # Find all positions where a record starts (with '>')
    record_starts = [0] + [m.start() for m in re.finditer(r'\n>', fasta_text)]

    # Extract each record
    records = []
    for i in range(len(record_starts)):
        start = record_starts[i]
        end = len(fasta_text) if i == len(record_starts) - 1 else record_starts[i + 1]
        record = fasta_text[start:end].strip()
        if record:
            records.append(record)

    return records


def run_hmmsearch_on_record(query_hmmfile: str, target_sequence: str, use_tempfile: bool = False, cpu: int = 1
                            ) -> List[SearchIO._model.query.QueryResult]:  # pylint: disable=protected-access
    """ Run hmmsearch on a HMM file and a single fasta record

        Arguments:
            query_hmmfile: the path to the HMM file
            target_sequence: a single fasta record as a string
            use_tempfile: if True, a tempfile will be written for the fasta input
                          instead of piping
            cpu: number of CPUs to use for this search (default: 1)

        Returns:
            a list of hmmsearch results as parsed by SearchIO
    """
    config = get_config()
    command = [config.executables.hmmsearch, "--cpu", str(cpu),
               "-o", os.devnull,  # throw away the verbose output
               "--domtblout", "result.domtab",
               query_hmmfile]

    with TemporaryDirectory(change=True):
        try:
            if use_tempfile:
                with open("input.fa", "w", encoding="utf-8") as handle:
                    handle.write(target_sequence)
                command.append("input.fa")
                run_result = execute(command)
            else:
                command.append('-')
                run_result = execute(command, stdin=target_sequence)
        except OSError:
            return []
        if not run_result.successful():
            logging.error('hmmsearch returned %d: %s while searching %s',
                          run_result.return_code, run_result.stderr, query_hmmfile)
            raise RuntimeError("Running hmmsearch failed.")
        return list(SearchIO.parse("result.domtab", 'hmmsearch3-domtab'))


def run_hmmsearch_parallel(query_hmmfile: str, target_sequence: str, use_tempfile: bool = False,
                           ) -> List[SearchIO._model.query.QueryResult]:  # pylint: disable=protected-access
    """ Run hmmsearch on a HMM file and a fasta input in parallel

        Splits the input fasta into individual records and processes each record
        in parallel, with each worker using only 1 CPU.

        Arguments:
            query_hmmfile: the path to the HMM file
            target_sequence: the fasta input to search as a string
            use_tempfile: if True, a tempfile will be written for the fasta input
                          instead of piping

        Returns:
            a list of hmmsearch results as parsed by SearchIO
    """
    if not target_sequence:
        return []

    # Split the fasta into individual records
    records = split_fasta(target_sequence)
    if not records:
        return []

    # Prepare arguments for parallel execution
    cpus = get_effective_cpus()
    cpu_per_process = 1
    args = [(query_hmmfile, record, use_tempfile, cpu_per_process) for record in records]

    # Run hmmsearch on each record in parallel
    results = parallel_function(run_hmmsearch_on_record, args, cpus=cpus)

    # Flatten the list of results
    flattened_results = []
    for result_list in results:
        flattened_results.extend(result_list)

    return flattened_results


def run_hmmsearch_version() -> str:
    """ Get the version of the hmmsearch """

    hmmsearch = get_config().executables.hmmsearch
    command = [
        hmmsearch,
        "-h",
    ]

    help_text = execute(command).stdout
    if not help_text.startswith("# hmmsearch"):
        msg = "unexpected output from hmmsearch: %s, check path"
        raise RuntimeError(msg % hmmsearch)

    version_line = help_text.split('\n')[1]
    return version_line.split()[2]
