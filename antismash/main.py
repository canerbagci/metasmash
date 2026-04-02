# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" The core setup functions wrapping all detection, analysis, and output
    modules together.

    The intended entry point of antismash is run_antismash() in this file.
"""

from collections import defaultdict
import cProfile
from datetime import datetime
import importlib
from io import StringIO
import glob
import logging
import os
import pkgutil
import pstats
import shutil
import time
import tempfile
import copy
import gc
import traceback
from typing import cast, Any, Dict, Iterator, List, Optional, Tuple, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from antismash.config import (
    ConfigType,
    get_config,
    update_config,
)
from antismash.common import errors, logs, memory as memory_diagnostics, record_processing, serialiser
from antismash.common.errors import AntismashInputError
from antismash.common.module_results import ModuleResults, DetectionResults
from antismash.common.path import get_full_path
from antismash.common.secmet import Record
from antismash.common import subprocessing
from antismash.detection import DetectionStage
from antismash.outputs import html
from antismash.support import genefinding
from antismash.custom_typing import AntismashModule

__version__ = "8.dev"


def _gather_analysis_modules() -> List[AntismashModule]:
    modules = []
    for module_data in pkgutil.walk_packages([get_full_path(__file__, "modules")]):
        module = importlib.import_module(f"antismash.modules.{module_data.name}")
        modules.append(cast(AntismashModule, module))
    return modules


def _gather_detection_modules() -> Dict[DetectionStage, List[AntismashModule]]:
    modules: Dict[DetectionStage, List[AntismashModule]] = {
        DetectionStage.FULL_GENOME: [],
        DetectionStage.AREA_FORMATION: [],
        DetectionStage.AREA_REFINEMENT: [],
        DetectionStage.PER_AREA: [],
    }
    for module_data in pkgutil.walk_packages([get_full_path(__file__, "detection")]):
        name = f"antismash.detection.{module_data.name}"
        module = cast(AntismashModule, importlib.import_module(name))
        stage = getattr(module, "DETECTION_STAGE", "")
        if not stage:
            raise ValueError(f"detection module missing DETECTION_STAGE attribute: {name}")
        assert isinstance(stage, DetectionStage)
        if stage not in modules:
            raise ValueError(f"detection module with unknown detection stage: {stage}")
        modules[stage].append(module)
    return modules


_ANALYSIS_MODULES = _gather_analysis_modules()
_DETECTION_MODULES = _gather_detection_modules()


def get_all_modules() -> List[AntismashModule]:
    """ Return a list of default modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    all_modules = []
    for modules in [get_detection_modules(), get_analysis_modules(),
                    get_output_modules(), get_support_modules()]:
        all_modules.extend(modules)
    return list(all_modules)


def get_detection_modules() -> List[AntismashModule]:
    """ Return a list of default detection modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    modules = []
    for stage in _DETECTION_MODULES.values():
        modules.extend(stage)
    return modules


def get_analysis_modules() -> List[AntismashModule]:
    """ Return a list of default analysis modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return list(_ANALYSIS_MODULES)


def get_output_modules() -> List[AntismashModule]:
    """ Return a list of default output modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    return [html]  # type: ignore  # a lot of casting avoided


def get_support_modules() -> List[AntismashModule]:
    """ Return a list of support modules

        Arguments:
            None

        Returns:
            a list of modules
    """
    genef = cast(AntismashModule, genefinding)
    return [genef]


def verify_options(options: ConfigType, modules: List[AntismashModule]) -> bool:
    """ Find and display any incompatibilities in provided options

        Arguments:
            options: the options to check
            modules: the modules to check the options of

        Returns:
            True if no problems detected, otherwise False
    """
    errors_found: List[str] = []
    configured_phase1_batch = int(getattr(options, "streaming_phase1_batch_size", 0) or 0)
    configured_phase2_window = int(getattr(options, "streaming_phase2_window_size", 0) or 0)
    if configured_phase1_batch < 0:
        errors_found.append("--streaming-phase1-batch-size must be >= 0")
    if configured_phase2_window < 0:
        errors_found.append("--streaming-phase2-window-size must be >= 0")
    for module in modules:
        try:
            logging.debug("Checking options for %s", module.__name__)
            errors_found.extend(module.check_options(options))
        except ValueError as err:
            errors_found.append(str(err))
    if not errors_found:
        return True

    logging.error("Incompatible options detected:\n  %s", "\n  ".join(errors_found))
    for error in errors_found:
        print(error)  # still commandline args, so don't use logging
    return False


def run_detection(record: Record, options: ConfigType,
                  module_results: Dict[str, Union[ModuleResults, Dict[str, Any]]]) -> Dict[str, float]:
    """ Detect different secondary metabolite clusters, PFAMs, and domains.

        Arguments:
            record: the Record to run detection over
            options: antiSMASH config
            module_results: a dictionary mapping a module's name to results from
                            a previous run on this module, as a ModuleResults subclass
                            or in JSON form

        Returns:
            the time taken by each detection module as a dictionary
    """
    timings: Dict[str, float] = {}

    try:
        # run full genome detections
        for module in _DETECTION_MODULES[DetectionStage.FULL_GENOME]:
            run_module(record, module, options, module_results, timings)
            results = module_results.get(module.__name__)
            if results:
                assert isinstance(results, ModuleResults)
                logging.debug("Adding detection results from %s to record %s", module.__name__, record.id)
                results.add_to_record(record)

        # generate cluster predictions
        logging.info("Detecting secondary metabolite clusters for record %s", record.id)
        modules = list(_DETECTION_MODULES[DetectionStage.AREA_FORMATION])
        modules.extend(_DETECTION_MODULES[DetectionStage.AREA_REFINEMENT])
        for module in modules:
            run_module(record, module, options, module_results, timings)
            results = module_results.get(module.__name__)
            if results:
                assert isinstance(results, DetectionResults), f"{module.__name__}, {type(results)}"
                for protocluster in results.get_predicted_protoclusters():
                    record.add_protocluster(protocluster)
                for region in results.get_predicted_subregions():
                    record.add_subregion(region)

        logging.debug("%d protoclusters found in record %s", len(record.get_protoclusters()), record.id)
        logging.debug("%d subregions found in record %s", len(record.get_subregions()), record.id)

        record.create_candidate_clusters()
        record.create_regions()

        if not record.get_regions():
            logging.info("No regions detected, skipping record %s", record.id)
            record.skip = "No regions detected"
            return timings

        logging.info("%d region(s) detected in record %s", len(record.get_regions()), record.id)

        # finally, run any detection limited to genes in clusters
        for module in _DETECTION_MODULES[DetectionStage.PER_AREA]:
            run_module(record, module, options, module_results, timings)
            results = module_results.get(module.__name__)
            if results:
                assert isinstance(results, ModuleResults)
                logging.debug("Adding detection results from %s to record %s", module.__name__, record.id)
                results.add_to_record(record)

        return timings
    except Exception as e:
        logging.error("Detection failed for record %s: %s", record.id, str(e))
        raise


def run_module(record: Record, module: AntismashModule, options: ConfigType,
               module_results: Dict[str, Union[ModuleResults, Dict[str, Any]]],
               timings: Dict[str, float]
               ) -> None:
    """ Run a module on a record

        Arguments:
            record: the record to run the analysis on
            module: the module to run, only run if enabled and not reusing results
            options: antismash Config
            module_results: a dictionary of module name to ModuleResults
                            instances or their JSON representations,
                            updated if the module runs
            timings: a dictionary mapping module name to time taken for that
                     module, will be updated with the module timing

        Returns:
            None
    """
    previous_results = module_results.pop(module.__name__, None)
    results = None
    if previous_results is not None:
        assert isinstance(previous_results, dict)
        logging.debug("Regenerating results for %s", module.__name__)
        results = module.regenerate_previous_results(previous_results, record, options)
        if results:
            module_results[module.__name__] = results
    # Check if all_enabled_modules exists in options (it won't in picklable options for workers)
    if hasattr(options, 'all_enabled_modules') and module not in options.all_enabled_modules:
        return
    assert results is None or isinstance(results, ModuleResults)

    logging.debug("Checking if %s should be run", module.__name__)
    if not module.is_enabled(options):
        return

    logging.info("Running %s on record %s", module.__name__, record.id)

    start = time.time()
    try:
        results = module.run_on_record(record, results, options)
    except Exception as e:
        logging.error("Error running %s on record %s: %s", module.__name__, record.id, str(e))
        raise
    duration = time.time() - start

    assert isinstance(results, ModuleResults), f"{module.__name__} returned {type(results)}"
    module_results[module.__name__] = results
    timings[module.__name__] = duration


def analyse_record(record: Record, options: ConfigType, modules: List[AntismashModule],
                   previous_result: Dict[str, Union[ModuleResults, Dict[str, Any]]]) -> Dict[str, float]:
    """ Run analysis modules on a record

        Arguments:
            record: the record to run the analysis on
            options: antismash Config
            modules: the modules to analyse with
                        each module will run only if enabled and not reusing all
                        results
            previous_result: a dictionary of module name to json results,
                                json results will be replaced by ModuleResults
                                instances

        Returns:
            a dictionary mapping module name to time taken
    """
    timings: Dict[str, float] = {}
    try:
        # try to run the given modules over the record
        for module in modules:
            run_module(record, module, options, previous_result, timings)
        return timings
    except Exception as e:
        logging.error("Analysis failed for record %s: %s", record.id, str(e))
        raise


def prepare_output_directory(name: str, input_file: str) -> None:
    """ Ensure the ouptut directory exists and is usable

        Raises an exception if the directory is unusable,
        or if results not being reused and directory not empty

        Arguments:
            name: the path of the directory
            input_file: the path of the input file

        Returns:
            None
    """
    # if not supplied, set the output directory to be the sequence name
    input_prefix = os.path.basename(canonical_base_filename(input_file, "", get_config()))
    if not name:
        name = os.path.abspath(input_prefix)
        update_config({"output_dir": name})

    if os.path.exists(name):
        if not os.path.isdir(name):
            raise AntismashInputError("Output directory {name!r} exists and is not a directory")
        # not empty (apart from a possible input dir), and not reusing its results
        if not input_file.endswith(".json") and \
                list(filter(_ignore_patterns, glob.glob(os.path.join(name, "*")))):
            raise AntismashInputError("Output directory contains other files, aborting for safety")

        # --reuse
        logging.debug("Removing existing region genbank files")
        for genbank in glob.glob(os.path.join(name, "*.region???.gbk")):
            os.remove(genbank)
        logging.debug("Reusing output directory: %s", name)
    else:
        logging.debug("Creating output directory: %s", name)
        os.mkdir(name)


def _ignore_patterns(entry: str) -> bool:
    """File name patterns that we want to ignore for the "outdir is empty" check."""
    config = get_config()
    if entry.endswith('/input') and os.path.isdir(entry):
        return False
    if os.path.abspath(entry) == os.path.abspath(config.logfile):
        return False

    return True


def write_profiling_results(profiler: cProfile.Profile, target: str) -> None:
    """ Write profiling files to file in human readable form and as a binary
        blob for external tool use (with the extra extension '.bin').

        If the file cannot be opened or written to, a shortened form will be
        written to stdout to avoid losing the data.

        Arguments:
            profiler: the profiler instance to log results of
            target: the path of the file to store reuslts in

        Returns:
            None
    """
    stream = StringIO()
    sortby = 'tottime'
    stats = pstats.Stats(profiler, stream=stream).sort_stats(sortby)
    stats.dump_stats(target + ".bin")
    stats.print_stats(.25)  # limit to the more meaningful first 25%
    stats.print_callers(.25)
    try:
        path_to_remove = os.path.dirname(os.path.realpath(__file__)) + os.path.sep
        with open(target, "w", encoding="utf-8") as handle:
            handle.write(stream.getvalue().replace(path_to_remove, ""))
        logging.info("Profiling report written to %s", target)
    except IOError:
        # if can't save to file, print to terminal, but only the head
        logging.debug("Couldn't open file to store profiling output")
        stream.truncate(0)
        stats.print_stats(20)  # first 20 lines only
        print(stream.getvalue())


def add_antismash_comments(records: List[Tuple[Record, SeqRecord]], options: ConfigType) -> None:
    """ Add antismash meta-annotation to records for genbank output

        Arguments:
            records: a list of Record, SeqRecord pairs
            options: antismash options

        Returns:
            None

    """
    if not records:
        return
    base_comment = {
        "Version": str(options.version),
        "Run date": str(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    }
    # include start/end details if relevant
    if options.start != -1 or options.end != -1:
        start = 1 if options.start == -1 else options.start
        # start/end is only valid for single records, as per record_processing
        if len(records) != 1:
            raise ValueError(
                "--start/--end options are only valid with a single input record, "
                f"but {len(records)} records were found"
            )
        end = len(records[0][0].seq) if options.end == -1 else options.end
        base_comment.update({
            "NOTE": "This is an extract from the original record!",
            "Starting at": str(start),
            "Ending at": str(end),
        })

    for record, bio_record in records:
        comment = dict(base_comment)
        if record.original_id:
            comment["Original ID"] = record.original_id

        if "structured_comment" not in bio_record.annotations:
            bio_record.annotations["structured_comment"] = {}

        bio_record.annotations["structured_comment"]["antiSMASH-Data"] = comment


def _filter_records_for_output(
        records: List,
        results: List[Dict[str, Any]],
        options: ConfigType,
) -> Tuple[List, List[Dict[str, Any]], int, List[str]]:
    """ Filter out records without regions when --skip-records-without-regions is set.

        Arguments:
            records: list of Record objects
            results: parallel list of per-record result dicts
            options: antismash Config instance

        Returns:
            (filtered_records, filtered_results, skipped_count, skipped_ids)
    """
    if not getattr(options, "output_skip_records_without_regions", False):
        return records, results, 0, []

    filtered_records = []
    filtered_results = []
    skipped_ids: List[str] = []

    for record, result in zip(records, results):
        if record.get_regions():
            filtered_records.append(record)
            filtered_results.append(result)
        else:
            skipped_ids.append(record.id)

    skipped_count = len(skipped_ids)
    if skipped_count:
        logging.info("Skipping %d records without regions from output: %s",
                     skipped_count, ", ".join(skipped_ids))

    return filtered_records, filtered_results, skipped_count, skipped_ids


def write_outputs(results: serialiser.AntismashResults, options: ConfigType) -> None:
    """ Write output files (webpage, genbank files, etc) to the output directory

        Arguments:
            results: a serialiser.AntismashResults instance
            options: an antismash Config instance

        Returns:
            None
    """
    logging.debug("Writing non-HTML output files")
    # don't use results for which the module no longer exists to regenerate/calculate
    module_results_per_record = []
    assert len(results.records) == len(results.results)
    for record, record_results in zip(results.records, results.results):
        record_result = {}
        for module_name, result in record_results.items():
            if isinstance(result, ModuleResults):
                assert result.record_id == record.id
                logging.debug(" Writing relevant output files for %s", module_name)
                record_result[module_name] = result
                result.write_outputs(record, options)
        module_results_per_record.append(record_result)

    # Apply --skip-records-without-regions filter
    filtered_records, filtered_module_results, skipped_count, _skipped_ids = \
        _filter_records_for_output(results.records, module_results_per_record, options)

    if html.is_enabled(options):
        logging.debug("Creating results page")
        start = time.time()
        html.write(filtered_records, filtered_module_results, options, get_all_modules(),
                   skipped_record_count=skipped_count)
        # use an average of times for html
        duration = (time.time() - start) / max(len(results.records), 1)
        for val in results.timings_by_record.values():
            val[html.__name__] = duration

    # convert records to biopython (only filtered records for output)
    bio_records = [record.to_biopython() for record in filtered_records]

    # add antismash meta-annotation to records
    add_antismash_comments(list(zip(filtered_records, bio_records)), options)

    if options.region_gbks:
        logging.debug("Writing cluster-specific genbank files")
        for record, bio_record in zip(filtered_records, bio_records):
            for region in record.get_regions():
                region.write_to_genbank(directory=options.output_dir, record=bio_record)

    # write records to an aggregate output
    base_filename = canonical_base_filename(results.input_file, options.output_dir, options)
    if options.summary_gbk:
        combined_filename = base_filename + ".gbk"
        logging.debug("Writing final genbank file to '%s'", combined_filename)
        SeqIO.write(bio_records, combined_filename, "genbank")

    zipfile = base_filename + ".zip"
    if os.path.exists(zipfile):
        os.remove(zipfile)
    if options.zip_output:
        logging.debug("Zipping output to '%s'", zipfile)
        with tempfile.NamedTemporaryFile(prefix="as_zip_tmp", suffix=".zip") as temp:
            shutil.make_archive(temp.name.replace(".zip", ""), "zip", root_dir=options.output_dir)
            shutil.copy(temp.name, zipfile)
            os.chmod(zipfile, 0o644)
        assert os.path.exists(zipfile)


def canonical_base_filename(input_file: str, directory: str, options: ConfigType) -> str:
    """Generate a canonical base filename if one isn't specified in the options."""
    if options.output_basename:
        base_filename = options.output_basename
    else:
        base_filename, ext = os.path.splitext(os.path.basename(input_file))
        if ext.lower() in (".gz", ".bz", ".xz"):
            base_filename, _ = os.path.splitext(base_filename)
        update_config({"output_basename": base_filename})

    return os.path.join(directory, base_filename)


def annotate_records(results: serialiser.AntismashResults) -> None:
    """ Annotates all analysed records with the results generated from them

        Arguments:
            results: a serialiser.AntismashResults instance

        Returns:
            None
    """
    detection_module_names = {mod.__name__ for mod in get_detection_modules()}
    for record, record_results in zip(results.records, results.results):
        if record.skip:
            logging.debug("Not annotating skipped record %s: %s", record.id, record.skip)
            continue
        if not record_results:
            logging.debug("No results for record %s, not annotating", record.id)
            continue
        logging.debug("Annotating record %s with results from: %s", record.id,
                      ", ".join([name.split()[0].split('.')[-1] for name in record_results]))
        for module, result in record_results.items():
            if module in detection_module_names:
                continue
            logging.debug(" Adding results from %s", module)
            assert isinstance(result, ModuleResults), type(result)
            result.add_to_record(record)


def read_data(sequence_file: Optional[str], options: ConfigType) -> serialiser.AntismashResults:
    """ Reads in the data to be used in the analysis run. Can be provided as
        as a sequence file (fasta/genbank) or as file of prior results

        Arguments:
            sequence_file: A fasta/genbank file to read (or None)
            options: An antismash Config instance

        Returns:
            a AntismashResults instance, populated only if reusing results

    """
    if not sequence_file and not options.reuse_results:
        raise ValueError("No sequence file or prior results to read")

    if sequence_file:
        records = record_processing.parse_input_sequence(
            sequence_file, options.taxon, options.minlength, options.start,
            options.end, gff_file=options.genefinding_gff3,
            ignore_invalid_records=not options.abort_on_invalid_records,
        )
        results = serialiser.AntismashResults(sequence_file.rsplit(os.sep, 1)[-1],
                                              records, [{} for i in records],
                                              __version__, taxon=options.taxon)
        update_config({"input_file": os.path.splitext(results.input_file)[1]})
    else:
        logging.debug("Attempting to reuse previous results in: %s", options.reuse_results)
        with open(options.reuse_results, "rb") as handle:
            contents = handle.read()
            if not contents:
                raise ValueError(f"No results contained in file: {options.reuse_results!r}")
        results = serialiser.AntismashResults.from_file(options.reuse_results)
        for record in results.records:
            record.strip_antismash_annotations()
        if options.taxon != results.taxon:
            logging.info("Reusing taxon %s from prior results", results.taxon)
            update_config({"taxon": results.taxon})

    update_config({"input_file": os.path.splitext(results.input_file)[0]})
    return results


def prepare_module_data(modules: Optional[List[AntismashModule]] = None) -> None:
    """ Calls all given modules' data preparation methods.

        Arguments:
            modules: a list of modules to use, if None all module will be used

        Returns:
            None
    """
    if modules is None:
        modules = get_all_modules()
    for module in modules:
        if hasattr(module, "prepare_data"):
            module.prepare_data()


def check_prerequisites(modules: List[AntismashModule], options: ConfigType) -> None:
    """ Checks that each module's prerequisites are satisfied. If not satisfied,
        a RuntimeError is raised.

        Any issues found are logged to logging.error

        Arguments:
            modules: the modules to check

        Returns:
            None
    """
    errors_by_module = {}
    for module in modules:
        logging.debug("Checking prerequisites for %s", module.__name__)
        res = module.check_prereqs(options)
        if res:
            errors_by_module[module.__name__] = res
    if errors_by_module:
        for module_name, errors_found in errors_by_module.items():
            for error in errors_found:
                logging.error("%s: preqrequisite failure: %s", module_name, error)
        raise RuntimeError("Modules failing prerequisites")


def list_plugins() -> None:
    """ Prints the name and short description of all registered modules

        Returns:
            None
    """

    def print_modules(modules: List[AntismashModule], indent: int = 0) -> None:
        max_length = max(len(mod.NAME) for mod in modules)
        format_string = f"{' ' * indent}%-{max_length}s:  %s"
        for module in modules:
            print(format_string % (module.NAME, module.SHORT_DESCRIPTION))
        print()

    print("Available plugins")
    print("  Detection modules")
    for stage, modules in _DETECTION_MODULES.items():
        simple_stage = str(stage).split(".")[1].replace("_", " ").capitalize()
        print(f"    {simple_stage}")
        print_modules(modules, indent=6)
    for title, modules in [
        ("Analysis", get_analysis_modules()),
        ("Output", get_output_modules()),
        ("Support", get_support_modules()),
    ]:
        print(f"  {title}")
        print_modules(modules, indent=4)


def log_module_runtimes(timings: Dict[str, Dict[str, float]]) -> None:
    """ Log the aggregate time taken per module.

        Arguments:
            timings: a dictionary mapping record id to
                        a dictionary mapping module name to time taken

        Returns:
            None
    """
    total_times: Dict[str, float] = defaultdict(lambda: 0.)
    for result in timings.values():
        for module, runtime in result.items():
            total_times[module] += runtime
    if not total_times:
        return
    logging.debug("Total times taken by modules")
    for module, runtime in sorted(total_times.items()):
        logging.debug("  %s: %.1fs", module, runtime)


def run_antismash(sequence_file: Optional[str], options: ConfigType) -> int:
    """ The complete antismash pipeline. Reads in data, runs detection and
        analysis modules over any records found, then outputs the results to
        file.

        Arguments:
            sequence_file: the sequence file to read in records from, can be
                            None if reusing results
            options: command line options

        Returns:
            0 if requested operations completed succesfully, otherwise 1
            Exceptions may also be raised
    """

    with logs.changed_logging(logfile=options.logfile, verbose=options.verbose,
                              debug=options.debug):
        try:
            result = _run_antismash(sequence_file, options)
        except errors.AntismashInputError as err:
            logging.error(str(err))
            raise
    return result


def _get_all_enabled_modules(modules: list[AntismashModule], options: ConfigType) -> list[AntismashModule]:
    return [module for module in modules if module.is_enabled(options)]


def process_record_detection(record_and_results: tuple, options: ConfigType) -> tuple:
    """ Process a single record with detection modules only.

        Designed to be called by parallel_function for record-level parallelism.

        Arguments:
            record_and_results: a tuple of (record, module_results)
            options: antiSMASH config

        Returns:
            a tuple of (record, module_results, timings)
    """
    record, module_results = record_and_results

    if record.skip:
        return record, module_results, {}

    logging.info("Starting detection for record: %s", record.id)
    try:
        timings = run_detection(record, options, module_results)
        logging.info("Completed detection for record: %s", record.id)
        return record, module_results, timings
    except Exception as e:
        logging.error("Error during detection for record %s: %s", record.id, str(e))
        raise


def process_record_analysis(record_and_results: tuple, options: ConfigType) -> tuple:
    """ Process a single record with analysis modules only.

        Designed to be called by parallel_function for record-level parallelism.

        Arguments:
            record_and_results: a tuple of (record, module_results)
            options: antiSMASH config

        Returns:
            a tuple of (record, module_results, timings)
    """
    record, module_results = record_and_results

    if record.skip or not record.get_regions():
        return record, module_results, {}

    logging.info("Starting analysis for record: %s", record.id)
    try:
        analysis_modules = get_analysis_modules()
        timings = analyse_record(record, options, analysis_modules, module_results)
        logging.info("Completed analysis for record: %s", record.id)
        return record, module_results, timings
    except Exception as e:
        logging.error("Error during analysis for record %s: %s", record.id, str(e))
        raise


_STREAMING_RECORD_THRESHOLD = 10


def should_use_streaming(options: ConfigType, sequence_file: Optional[str],
                         ) -> Tuple[bool, Optional[List[Tuple[str, int]]]]:
    """ Determine whether to use the streaming pipeline.

        Arguments:
            options: antismash config
            sequence_file: the input sequence file path

        Returns:
            a tuple of (use_streaming, prefetched_metadata) where
            prefetched_metadata is a list of (record_id, seq_length) tuples
            when auto mode scanned them (to avoid re-scanning), or None
    """
    streaming = getattr(options, 'streaming', 'auto')

    if streaming == 'off':
        return False, None
    if streaming == 'on':
        if options.reuse_results:
            logging.warning("Streaming mode cannot be used with --reuse-results, falling back to classic")
            return False, None
        if options.genefinding_gff3:
            logging.warning("Streaming mode cannot be used with --genefinding-gff3, falling back to classic")
            return False, None
        return True, None

    # auto mode
    if not sequence_file or options.reuse_results:
        return False, None
    if options.genefinding_gff3:
        return False, None

    # scan record metadata (lightweight — no sequence data retained)
    try:
        metadata = record_processing.scan_record_metadata(
            sequence_file, options.minlength,
            ignore_invalid_records=not getattr(options, 'abort_on_invalid_records', True),
        )
    except Exception:
        return False, None  # if scanning fails, fall back to classic

    if len(metadata) > _STREAMING_RECORD_THRESHOLD:
        logging.info("Auto-enabling streaming mode for %d records (threshold: %d)",
                     len(metadata), _STREAMING_RECORD_THRESHOLD)
        return True, metadata
    return False, None


def process_record_full(record_tuple: tuple, options: ConfigType) -> tuple:
    """ Full pipeline for a single record. Runs in a parallel worker.

        Combines secmet conversion, sanitisation, gene finding, detection,
        and analysis into a single worker function.

        Arguments:
            record_tuple: a tuple of (SeqRecord, record_index)
            options: antismash config

        Returns:
            a tuple of (Record, module_results_dict, timings_dict)
    """
    import functools

    bio_record, record_index = record_tuple

    # 1. Strip old antiSMASH annotations
    record_processing.strip_record(bio_record)

    # 2. Convert to secmet Record
    try:
        record = Record.from_biopython(bio_record, options.taxon, discard_antismash_features=True)
    except Exception as err:
        logging.error("Failed to convert record %s to secmet: %s", bio_record.id, err)
        raise

    record.strip_antismash_annotations()
    record.record_index = record_index

    # preserve original_id if set on the SeqRecord
    original_id = getattr(bio_record, 'original_id', None)
    if original_id and not record.original_id:
        record.original_id = original_id

    # 3. Sanitise sequence
    record = record_processing.sanitise_sequence(record)
    if record.skip:
        return record, {}, {}

    # check for empty records
    if not record.seq:
        logging.warning("Record %s has no sequence, skipping.", record.id)
        record.skip = "contains no sequence"
        return record, {}, {}

    if not record.id:
        logging.error("Record has no name")
        record.skip = "no name"
        return record, {}, {}

    # 4. Gene finding
    genefinding_opts = {key: val for key, val in options if key.startswith("genefinding")}
    genefinding_opts["taxon"] = options.taxon
    base_function = genefinding.run_on_record
    partial = functools.partial(record_processing.ensure_cds_info, base_function,
                                **genefinding_opts)
    record = partial(record)
    if record.skip:
        return record, {}, {}

    # 5. Detection
    module_results: Dict[str, Any] = {}
    try:
        timings = run_detection(record, options, module_results)
    except Exception as e:
        logging.error("Detection failed for record %s: %s", record.id, str(e))
        raise

    # 6. Analysis (if regions found)
    if record.get_regions():
        try:
            analysis_timings = analyse_record(record, options,
                                              get_analysis_modules(), module_results)
            timings.update(analysis_timings)
        except Exception as e:
            logging.error("Analysis failed for record %s: %s", record.id, str(e))
            raise

    return record, module_results, timings


_RECORD_FAILED = "__RECORD_PROCESSING_FAILED__"  # sentinel for records that failed in parallel processing


def _format_worker_memory_delta(before_mb: Optional[float], after_mb: Optional[float]) -> str:
    """ Format an RSS delta for diagnostics logging. """
    if before_mb is None or after_mb is None:
        return "n/a"
    return f"{after_mb - before_mb:+.1f} MB"


def _log_worker_memory(label: str, record: Optional[Record], options: ConfigType,
                       rss_before_mb: Optional[float],
                       record_id: Optional[str] = None,
                       timings: Optional[Dict[str, float]] = None,
                       error: Optional[str] = None) -> None:
    """ Log worker-side memory diagnostics for a single record. """
    if not memory_diagnostics.diagnostics_enabled(options):
        return
    resolved_record_id = record.id if record is not None else record_id or "unknown"
    regions = len(record.get_regions()) if record is not None else -1
    cds_count = len(record.get_cds_features()) if record is not None else -1
    total_time = sum(timings.values()) if timings else 0.0
    rss_after_mb = memory_diagnostics.current_rss_mb()
    logging.info(
        "memdiag worker phase=%s pid=%d record=%s rss=%s peak=%s delta=%s regions=%d cds=%d time=%.2fs%s",
        label,
        os.getpid(),
        resolved_record_id,
        memory_diagnostics.format_mb(rss_after_mb),
        memory_diagnostics.format_mb(memory_diagnostics.peak_rss_mb()),
        _format_worker_memory_delta(rss_before_mb, rss_after_mb),
        regions,
        cds_count,
        total_time,
        f" error={error}" if error else "",
    )


def _log_parent_memory(label: str, options: ConfigType,
                       extra: Optional[Dict[str, Union[int, str]]] = None,
                       trace_snapshot: Any = None) -> Any:
    """ Log parent-side memory diagnostics for streaming runs. """
    if not memory_diagnostics.diagnostics_enabled(options):
        return trace_snapshot

    stats = memory_diagnostics.process_tree_stats()
    counts = memory_diagnostics.tracked_object_counts()
    details = [
        f"parent_rss={memory_diagnostics.format_mb(stats['parent_rss_mb'])}",
        f"descendants={int(stats['descendant_count'] or 0)}",
        f"descendant_rss={memory_diagnostics.format_mb(stats['descendant_rss_mb'])}",
        f"known_descendants={int(stats['descendant_rss_known_count'] or 0)}",
        f"records={counts['Record']}",
        f"record_layers={counts['RecordLayer']}",
        f"region_layers={counts['RegionLayer']}",
        f"module_results={counts['ModuleResults']}",
    ]
    if extra:
        details.extend(f"{key}={value}" for key, value in extra.items())
    logging.info("memdiag parent %s %s", label, " ".join(details))
    if memory_diagnostics.tracemalloc_enabled(options):
        return memory_diagnostics.maybe_log_tracemalloc(label, trace_snapshot)
    return trace_snapshot


def _safe_process_record_full(record_tuple, options):
    """ Wrapper around process_record_full that catches exceptions for
        individual records, allowing the rest of the pool to continue.

        When abort_on_invalid_records is True, exceptions propagate normally.
        Otherwise, failures are logged and a sentinel tuple is returned.
    """
    try:
        return process_record_full(record_tuple, options)
    except Exception as e:
        bio_record, record_index = record_tuple
        if getattr(options, 'abort_on_invalid_records', True):
            raise
        tb = traceback.format_exc()
        logging.error("Record %s (index %d) failed, skipping: %s",
                      bio_record.id, record_index, e)
        return _RECORD_FAILED, bio_record.id, tb


def process_record_detection_streaming(record_tuple: tuple, options: ConfigType) -> tuple:
    """ Detection-only pipeline for a single record in streaming mode.

        Like process_record_full but stops after detection — no analysis step.
        Used in Phase 1 of two-phase streaming.

        Arguments:
            record_tuple: a tuple of (SeqRecord, record_index)
            options: antismash config

        Returns:
            a tuple of (Record, module_results_dict, timings_dict)
    """
    import functools

    bio_record, record_index = record_tuple

    # 1. Strip old antiSMASH annotations
    record_processing.strip_record(bio_record)

    # 2. Convert to secmet Record
    try:
        record = Record.from_biopython(bio_record, options.taxon, discard_antismash_features=True)
    except Exception as err:
        logging.error("Failed to convert record %s to secmet: %s", bio_record.id, err)
        raise

    record.strip_antismash_annotations()
    record.record_index = record_index

    # preserve original_id if set on the SeqRecord
    original_id = getattr(bio_record, 'original_id', None)
    if original_id and not record.original_id:
        record.original_id = original_id

    # 3. Sanitise sequence
    record = record_processing.sanitise_sequence(record)
    if record.skip:
        return record, {}, {}

    if not record.seq:
        logging.warning("Record %s has no sequence, skipping.", record.id)
        record.skip = "contains no sequence"
        return record, {}, {}

    if not record.id:
        logging.error("Record has no name")
        record.skip = "no name"
        return record, {}, {}

    # 4. Gene finding
    genefinding_opts = {key: val for key, val in options if key.startswith("genefinding")}
    genefinding_opts["taxon"] = options.taxon
    base_function = genefinding.run_on_record
    partial = functools.partial(record_processing.ensure_cds_info, base_function,
                                **genefinding_opts)
    record = partial(record)
    if record.skip:
        return record, {}, {}

    # 5. Detection (NO analysis)
    module_results: Dict[str, Any] = {}
    try:
        timings = run_detection(record, options, module_results)
    except Exception as e:
        logging.error("Detection failed for record %s: %s", record.id, str(e))
        raise

    return record, module_results, timings


def _safe_process_record_detection_streaming(record_tuple, options):
    """ Safe wrapper for process_record_detection_streaming.
        Returns _RECORD_FAILED sentinel on failure when abort_on_invalid_records is False.
    """
    rss_before_mb = memory_diagnostics.current_rss_mb() if memory_diagnostics.diagnostics_enabled(options) else None
    try:
        result = process_record_detection_streaming(record_tuple, options)
        record, _, timings = result
        _log_worker_memory("detection", record, options, rss_before_mb, timings=timings)
        return result
    except Exception as e:
        bio_record, record_index = record_tuple
        if getattr(options, 'abort_on_invalid_records', True):
            raise
        tb = traceback.format_exc()
        logging.error("Record %s (index %d) failed during detection, skipping: %s",
                      bio_record.id, record_index, e)
        _log_worker_memory("detection-failed", None, options, rss_before_mb,
                           record_id=bio_record.id, error=str(e))
        return _RECORD_FAILED, bio_record.id, tb


def _safe_process_record_analysis(record_and_results, options):
    """ Safe wrapper for process_record_analysis.
        Returns _RECORD_FAILED sentinel on failure when abort_on_invalid_records is False.
    """
    rss_before_mb = memory_diagnostics.current_rss_mb() if memory_diagnostics.diagnostics_enabled(options) else None
    try:
        result = process_record_analysis(record_and_results, options)
        record, _, timings = result
        _log_worker_memory("analysis", record, options, rss_before_mb, timings=timings)
        return result
    except Exception as e:
        record, _ = record_and_results
        if getattr(options, 'abort_on_invalid_records', True):
            raise
        tb = traceback.format_exc()
        logging.error("Record %s failed during analysis, skipping: %s",
                      record.id, e)
        _log_worker_memory("analysis-failed", record, options, rss_before_mb,
                           record_id=record.id, error=str(e))
        return _RECORD_FAILED, record.id, tb


def _strip_record_for_overview(record: Record) -> None:
    """Strip heavy data from a Record after all per-record output is written.

    Keeps record metadata (id, name, record_index, annotations), regions,
    protoclusters, candidate clusters, subregions, and CDS features (referenced
    by region.cds_children) so that RecordLayer/RegionLayer can still render
    overview.html and the dashboard.

    Removes: DNA sequence, CDS translations, and record-level feature indexes
    that are no longer needed.
    """
    # Cache GC content before clearing sequence (used by serialiser)
    try:
        record.get_gc_content()
    except (ValueError, ZeroDivisionError):
        pass  # empty or undefined sequence, nothing to cache

    # Clear the DNA sequence — typically the single largest item (3-6 MB)
    record.seq = Seq("")

    # Clear CDS translations (protein sequences, 1-2 MB total)
    for cds in record.get_cds_features():
        cds._translation = ""  # pylint: disable=protected-access

    # Clear record-level feature indexes (frees objects not referenced by regions)
    record._nonspecific_features.clear()  # pylint: disable=protected-access
    record._genes.clear()  # pylint: disable=protected-access
    record._genes_by_name.clear()  # pylint: disable=protected-access
    record._cds_by_name.clear()  # pylint: disable=protected-access
    record._cds_by_location.clear()  # pylint: disable=protected-access
    record._pfam_domains.clear()  # pylint: disable=protected-access
    record._pfams_by_cds_name.clear()  # pylint: disable=protected-access
    record._antismash_domains.clear()  # pylint: disable=protected-access
    record._antismash_domains_by_tool.clear()  # pylint: disable=protected-access
    record._antismash_domains_by_cds_name.clear()  # pylint: disable=protected-access
    record._domains_by_name.clear()  # pylint: disable=protected-access
    record._cds_motifs.clear()  # pylint: disable=protected-access
    record._modules.clear()  # pylint: disable=protected-access


def _preload_pfam_caches(options: ConfigType) -> None:
    """ Pre-load PFAM databases into module-level caches for fork CoW sharing.

        Arguments:
            options: antismash config
    """
    from antismash.common import pfamdb

    pfamdb.init_shared_cache()

    if hasattr(options, 'clusterhmmer') and options.clusterhmmer:
        database_version = options.clusterhmmer_pfamdb_version
        if database_version == "latest":
            database_version = pfamdb.find_latest_database_version(options.database_dir)
        database = os.path.join(options.database_dir, 'pfam', database_version, 'Pfam-A.hmm')
        pfamdb.preload_pfam_cutoffs(database)
        pfamdb.preload_pfam_mappings(database)

    if hasattr(options, 'fullhmmer') and options.fullhmmer:
        database_version = options.fullhmmer_pfamdb_version
        if database_version == "latest":
            database_version = pfamdb.find_latest_database_version(options.database_dir)
        database = os.path.join(options.database_dir, 'pfam', database_version, 'Pfam-A.hmm')
        pfamdb.preload_pfam_cutoffs(database)
        pfamdb.preload_pfam_mappings(database)


def _preload_analysis_caches(options: ConfigType) -> None:
    """ Pre-load analysis module databases into module-level caches for fork CoW sharing.

        Must be called before forking worker processes so that workers inherit
        the populated caches via copy-on-write, avoiding per-worker duplication.

        Arguments:
            options: antismash config
    """
    if hasattr(options, 'cb_general') and (
            options.cb_general or options.cb_subclusters or options.cb_knownclusters):
        logging.info("Pre-loading clusterblast databases for fork CoW sharing")
        from antismash.modules.clusterblast.core import preload_clusterblast_databases
        preload_clusterblast_databases(options)

    if hasattr(options, 'cc_mibig') and (options.cc_mibig or options.cc_custom_dbs):
        logging.info("Pre-loading cluster_compare databases for fork CoW sharing")
        from antismash.modules.cluster_compare import preload_databases
        preload_databases(options)

    if hasattr(options, 'pfam2go') and options.pfam2go:
        logging.info("Pre-loading pfam2go mapping for fork CoW sharing")
        from antismash.modules.pfam2go.pfam2go import preload_pfam2go_mapping
        preload_pfam2go_mapping()

    if hasattr(options, 'tfbs') and options.tfbs:
        logging.info("Pre-loading TFBS matrices for fork CoW sharing")
        from antismash.modules.tfbs_finder.tfbs_finder import preload_matrices
        preload_matrices()


_STREAMING_PHASE2_WINDOW_MIN = 1024


def _compute_streaming_phase1_batch_size(options: ConfigType) -> int:
    """Return the bounded Phase 1 detection batch size for streaming runs."""
    configured = int(getattr(options, "streaming_phase1_batch_size", 0) or 0)
    if configured < 0:
        raise ValueError("streaming_phase1_batch_size must be >= 0")
    if configured:
        return configured
    return max(_STREAMING_PHASE2_WINDOW_MIN, int(options.workers) * 4)


def _compute_streaming_phase2_window_size(options: ConfigType) -> int:
    """Return the bounded Phase 2 window size for streaming runs."""
    configured = int(getattr(options, "streaming_phase2_window_size", 0) or 0)
    if configured < 0:
        raise ValueError("streaming_phase2_window_size must be >= 0")
    if configured:
        return configured
    return max(_STREAMING_PHASE2_WINDOW_MIN, int(options.workers) * 4)


def _take_record_batch(record_iterator: Iterator[Tuple[SeqRecord, int]],
                       batch_size: int) -> List[Tuple[SeqRecord, int]]:
    """Take up to batch_size records from an iterator."""
    batch: List[Tuple[SeqRecord, int]] = []
    for _ in range(batch_size):
        try:
            batch.append(next(record_iterator))
        except StopIteration:
            break
    return batch


def _run_phase2_window(phase2_inputs: Dict[str, Tuple[Record, Dict[str, Any]]],
                       options: ConfigType,
                       picklable_options_p2: ConfigType,
                       user_workers: int,
                       json_writer: serialiser.StreamingJsonWriter,
                       all_modules: List[AntismashModule],
                       options_layer: Any,
                       data_writer: Any,
                       gbk_handle: Any,
                       timings_by_record: Dict[str, Dict[str, float]],
                       lightweight_records: List[Dict[str, Any]],
                       record_summaries: List[Any],
                       detection_module_names: set[str],
                       phase2_seen: int,
                       regions_count: int,
                       failed_count: int,
                       trace_snapshot: Any,
                       window_index: int,
                       window_size: int) -> Tuple[int, int, int, Any]:
    """Run one bounded Phase 2 analysis window."""
    from antismash.common.subprocessing import parallel_function_lazy
    from antismash.outputs.html.generator import generate_region_files_for_record

    effective_threads = max(1, options.cpus // user_workers)
    logging.info("Phase 2 window %d: Analysis of %d records with %d workers x %d threads",
                 window_index, len(phase2_inputs), user_workers, effective_threads)
    trace_snapshot = _log_parent_memory(
        "phase2-window-start",
        options,
        extra={
            "window_index": window_index,
            "window_size": window_size,
            "pending_phase2": len(phase2_inputs),
            "workers": user_workers,
            "threads_per_worker": effective_threads,
            "stored_record_summaries": len(record_summaries),
            "stored_lightweight_records": len(lightweight_records),
        },
        trace_snapshot=trace_snapshot,
    )

    def _analysis_args() -> Iterator[Tuple[Tuple[Record, Dict[str, Any]], ConfigType]]:
        record_ids = list(phase2_inputs.keys())
        for rec_id in record_ids:
            if rec_id in phase2_inputs:
                rec, mr = phase2_inputs[rec_id]
                yield ((rec, mr), picklable_options_p2)

    for item in parallel_function_lazy(
            _safe_process_record_analysis, _analysis_args(),
            cpus=user_workers):
        phase2_seen += 1
        if item[0] == _RECORD_FAILED:
            _, record_id, error_msg = item
            logging.warning("Analysis failed for record %s (see debug log)", record_id)
            logging.debug("Traceback for failed record %s:\n%s", record_id, error_msg)
            if record_id in phase2_inputs:
                rec, mr = phase2_inputs.pop(record_id)
                json_writer.write_record(rec, mr)
            failed_count += 1
            if (memory_diagnostics.diagnostics_enabled(options)
                    and phase2_seen % memory_diagnostics.diagnostics_interval(options) == 0):
                trace_snapshot = _log_parent_memory(
                    "phase2-progress",
                    options,
                    extra={
                        "window_index": window_index,
                        "window_size": window_size,
                        "seen": phase2_seen,
                        "completed": regions_count,
                        "pending_phase2": len(phase2_inputs),
                        "stored_record_summaries": len(record_summaries),
                        "stored_lightweight_records": len(lightweight_records),
                        "failed": failed_count,
                    },
                    trace_snapshot=trace_snapshot,
                )
            continue

        record, mod_results, analysis_timings = item

        if analysis_timings and record.id in timings_by_record:
            timings_by_record[record.id].update(analysis_timings)
        elif analysis_timings:
            timings_by_record[record.id] = analysis_timings

        filtered_results: Dict[str, ModuleResults] = {}
        for module_name, result in mod_results.items():
            if isinstance(result, ModuleResults):
                filtered_results[module_name] = result
                if module_name not in detection_module_names:
                    result.write_outputs(record, options)

        for module_name, result in filtered_results.items():
            if module_name not in detection_module_names:
                result.add_to_record(record)

        json_writer.write_record(record, mod_results)

        if options.region_gbks or gbk_handle is not None:
            bio_record = record.to_biopython()
            add_antismash_comments([(record, bio_record)], options)
            if options.region_gbks:
                for region in record.get_regions():
                    region.write_to_genbank(directory=options.output_dir,
                                            record=bio_record)
            if gbk_handle is not None:
                SeqIO.write([bio_record], gbk_handle, "genbank")
            del bio_record

        if options_layer is not None:
            light_record, record_summary = generate_region_files_for_record(
                record, filtered_results, options, all_modules,
                options_layer, data_writer=data_writer,
            )
            if record_summary is not None:
                lightweight_records.append(light_record)
                record_summaries.append(record_summary)

        regions_count += 1
        _strip_record_for_overview(record)
        phase2_inputs.pop(record.id, None)

        if (memory_diagnostics.diagnostics_enabled(options)
                and phase2_seen % memory_diagnostics.diagnostics_interval(options) == 0):
            trace_snapshot = _log_parent_memory(
                "phase2-progress",
                options,
                extra={
                    "window_index": window_index,
                    "window_size": window_size,
                    "seen": phase2_seen,
                    "completed": regions_count,
                    "pending_phase2": len(phase2_inputs),
                    "stored_record_summaries": len(record_summaries),
                    "stored_lightweight_records": len(lightweight_records),
                    "failed": failed_count,
                },
                trace_snapshot=trace_snapshot,
            )

    trace_snapshot = _log_parent_memory(
        "phase2-window-complete",
        options,
        extra={
            "window_index": window_index,
            "window_size": window_size,
            "seen": phase2_seen,
            "completed": regions_count,
            "pending_phase2": len(phase2_inputs),
            "stored_record_summaries": len(record_summaries),
            "stored_lightweight_records": len(lightweight_records),
            "failed": failed_count,
        },
        trace_snapshot=trace_snapshot,
    )
    return phase2_seen, regions_count, failed_count, trace_snapshot


def _run_antismash_streaming(sequence_file: str, options: ConfigType,
                             prefetched_metadata: Optional[List[Tuple[str, int]]] = None,
                             ) -> int:
    """ Memory-bounded streaming pipeline.

        Two-pass input: first pass collects lightweight (id, length) metadata,
        second pass streams only accepted records one at a time.  Parallel
        processing uses lazy iteration so only one result is in memory at a
        time on the consumer side.  Records without BGC regions are serialized
        to JSON immediately and discarded; only the ~0.4 % of records with
        regions are kept for HTML generation and GenBank output.

        Arguments:
            sequence_file: path to the input sequence file
            options: antismash config
            prefetched_metadata: optional metadata from should_use_streaming
                                 auto-mode (avoids re-scanning)

        Returns:
            0 on success, 1 if all records failed
    """
    from antismash.common.subprocessing import parallel_function_lazy
    from antismash.common.layers import OptionsLayer
    from antismash.outputs.html.generator import (
        finalize_streaming_html_output,
        StreamingRegionDataWriter,
    )
    from antismash.outputs.html.taxonomy import parse_taxonomy_file

    start_time = datetime.now()

    # --- Pass 1: scan metadata (id, length) without retaining sequences ---
    if prefetched_metadata is not None:
        metadata = prefetched_metadata
    else:
        metadata = record_processing.scan_record_metadata(
            sequence_file, options.minlength,
            ignore_invalid_records=not options.abort_on_invalid_records,
        )
    update_config({"input_file": os.path.splitext(os.path.basename(sequence_file))[0]})

    # Resolve IDs (dedup, shorten, apply --limit_to_record / --limit)
    accepted_ids, accepted_count = record_processing.resolve_record_ids(metadata, options)
    del metadata  # free ~660 MB for 11M records

    if accepted_count == 0:
        raise AntismashInputError("no records remaining after filtering")

    # Prepare output directory and HTML assets
    prepare_output_directory(options.output_dir, sequence_file)

    html_enabled = html.is_enabled(options)
    if html_enabled:
        html.copy_template_dir('css', options.output_dir, pattern=f"{options.taxon}.css")
        html.copy_template_dir('js', options.output_dir)
        local_js = os.path.join(options.output_dir, "js", "antismash.js")
        if not os.path.exists(local_js):
            js_path = html.find_local_antismash_js_path(options)
            if js_path:
                logging.debug("Results page using antismash.js from local copy: %s", js_path)
                shutil.copy(js_path, local_js)
        if not os.path.exists(local_js):
            logging.debug("Results page using antismash.js from remote host")
        html.copy_template_dir('images', options.output_dir)

    # Save user's --workers for Phase 2 (analysis); Phase 1 uses cpus workers
    user_workers = options.workers

    # --- Phase 1 setup: preload only PFAM caches, configure for many single-threaded workers ---
    update_config({"workers": options.cpus})  # get_effective_cpus() will return 1 in workers
    picklable_options = _create_picklable_options(options)
    _preload_pfam_caches(options)
    # Do NOT preload analysis caches yet — keeps Phase 1 workers lightweight
    gc.freeze()

    trace_snapshot = None
    memory_diagnostics.maybe_start_tracemalloc(options)
    phase1_batch_size = _compute_streaming_phase1_batch_size(options)
    phase2_window_size = _compute_streaming_phase2_window_size(options)

    # --- Pass 2: two-phase parallel processing ---
    logging.info("Processing %d records in two-phase streaming mode", accepted_count)
    accepted_records = iter(record_processing.iter_accepted_records(sequence_file, accepted_ids))

    # Open the streaming JSON writer
    input_basename = os.path.basename(sequence_file)
    json_filename = canonical_base_filename(input_basename, options.output_dir, options) + ".json"

    all_modules = get_all_modules()
    options_layer = OptionsLayer(options, all_modules) if html_enabled else None

    # Accumulators
    timings_by_record: Dict[str, Dict[str, float]] = {}

    # Lightweight accumulators for finalize_streaming_html_output
    lightweight_records: List[Dict[str, Any]] = []
    record_summaries: List[Any] = []

    total_processed = 0
    total_without_regions = 0
    regions_count = 0
    failed_count = 0
    skip_without_regions = getattr(options, "output_skip_records_without_regions", False)
    phase1_seen = 0
    phase2_seen = 0

    # Collect records with regions for Phase 2 (dict for O(1) error-recovery lookup)
    phase2_inputs: Dict[str, Tuple[Record, Dict[str, Any]]] = {}
    detection_module_names = {m.__name__ for m in get_detection_modules()}
    phase1_batch_index = 0
    window_index = 0
    source_exhausted = False
    analysis_caches_preloaded = False
    picklable_options_p2: Optional[ConfigType] = None
    data_handle = None
    data_writer = None
    gbk_handle = None

    logging.debug("Writing streaming json results to '%s'", json_filename)
    trace_snapshot = _log_parent_memory(
        "streaming-start",
        options,
        extra={
            "accepted_count": accepted_count,
            "workers": options.cpus,
            "analysis_workers": user_workers,
            "phase1_batch_size": phase1_batch_size,
            "window_size": phase2_window_size,
        },
        trace_snapshot=trace_snapshot,
    )
    with open(json_filename, "w", encoding="utf-8") as json_handle:
        json_writer = serialiser.StreamingJsonWriter(
            json_handle, input_basename, __version__, options.taxon)
        try:
            while not source_exhausted or phase2_inputs:
                while not source_exhausted and len(phase2_inputs) < phase2_window_size:
                    record_batch = _take_record_batch(accepted_records, phase1_batch_size)
                    if not record_batch:
                        source_exhausted = True
                        break

                    phase1_batch_index += 1
                    update_config({"workers": options.cpus})
                    logging.info("Phase 1 batch %d: Detection with %d workers x 1 thread over %d records",
                                 phase1_batch_index, options.cpus, len(record_batch))
                    trace_snapshot = _log_parent_memory(
                        "phase1-batch-start",
                        options,
                        extra={
                            "batch_index": phase1_batch_index,
                            "batch_size": len(record_batch),
                            "phase1_batch_size": phase1_batch_size,
                            "window_size": phase2_window_size,
                            "pending_phase2": len(phase2_inputs),
                            "processed": total_processed,
                            "without_regions": total_without_regions,
                            "failed": failed_count,
                        },
                        trace_snapshot=trace_snapshot,
                    )

                    def _record_args(batch: List[Tuple[SeqRecord, int]] = record_batch
                                     ) -> Iterator[Tuple[Tuple[SeqRecord, int], ConfigType]]:
                        for rec_tuple in batch:
                            yield (rec_tuple, picklable_options)

                    for item in parallel_function_lazy(
                            _safe_process_record_detection_streaming, _record_args(),
                            cpus=options.cpus):
                        phase1_seen += 1

                        if item[0] == _RECORD_FAILED:
                            _, record_id, error_msg = item
                            logging.warning("Skipping failed record %s (see debug log)", record_id)
                            logging.debug("Traceback for failed record %s:\n%s", record_id, error_msg)
                            failed_count += 1
                            if (memory_diagnostics.diagnostics_enabled(options)
                                    and phase1_seen % memory_diagnostics.diagnostics_interval(options) == 0):
                                trace_snapshot = _log_parent_memory(
                                    "phase1-progress",
                                    options,
                                    extra={
                                        "batch_index": phase1_batch_index,
                                        "phase1_batch_size": phase1_batch_size,
                                        "window_size": phase2_window_size,
                                        "seen": phase1_seen,
                                        "processed": total_processed,
                                        "pending_phase2": len(phase2_inputs),
                                        "without_regions": total_without_regions,
                                        "failed": failed_count,
                                    },
                                    trace_snapshot=trace_snapshot,
                                )
                            continue

                        record, mod_results, rec_timings = item
                        total_processed += 1
                        if rec_timings:
                            timings_by_record[record.id] = rec_timings

                        for module_name, result in mod_results.items():
                            if isinstance(result, ModuleResults):
                                assert result.record_id == record.id
                                result.write_outputs(record, options)

                        if record.get_regions():
                            phase2_inputs[record.id] = (record, mod_results)
                        elif not skip_without_regions:
                            json_writer.write_record(record, mod_results)
                            total_without_regions += 1
                        else:
                            total_without_regions += 1

                        if (memory_diagnostics.diagnostics_enabled(options)
                                and phase1_seen % memory_diagnostics.diagnostics_interval(options) == 0):
                            trace_snapshot = _log_parent_memory(
                                "phase1-progress",
                                options,
                                extra={
                                    "batch_index": phase1_batch_index,
                                    "phase1_batch_size": phase1_batch_size,
                                    "window_size": phase2_window_size,
                                    "seen": phase1_seen,
                                    "processed": total_processed,
                                    "pending_phase2": len(phase2_inputs),
                                    "without_regions": total_without_regions,
                                    "failed": failed_count,
                                },
                                trace_snapshot=trace_snapshot,
                            )

                    trace_snapshot = _log_parent_memory(
                        "phase1-batch-complete",
                        options,
                        extra={
                            "batch_index": phase1_batch_index,
                            "batch_size": len(record_batch),
                            "phase1_batch_size": phase1_batch_size,
                            "window_size": phase2_window_size,
                            "seen": phase1_seen,
                            "processed": total_processed,
                            "pending_phase2": len(phase2_inputs),
                            "without_regions": total_without_regions,
                            "failed": failed_count,
                        },
                        trace_snapshot=trace_snapshot,
                    )
                    gc.collect()

                    if len(phase2_inputs) >= phase2_window_size:
                        break

                if phase2_inputs:
                    window_index += 1
                    update_config({"workers": user_workers})
                    if picklable_options_p2 is None:
                        picklable_options_p2 = _create_picklable_options(options)
                    if not analysis_caches_preloaded:
                        _preload_analysis_caches(options)
                        gc.freeze()
                        analysis_caches_preloaded = True

                    if html_enabled and data_writer is None:
                        data_file_path = os.path.join(options.output_dir, "regions_data.js")
                        data_handle = open(data_file_path, "w", encoding="utf-8")
                        data_writer = StreamingRegionDataWriter(data_handle)

                    if options.summary_gbk and gbk_handle is None:
                        combined_filename = canonical_base_filename(
                            input_basename, options.output_dir, options
                        ) + ".gbk"
                        logging.debug("Writing final genbank file to '%s'", combined_filename)
                        gbk_handle = open(combined_filename, "w")

                    if options.region_gbks and window_index == 1:
                        logging.debug("Writing cluster-specific genbank files")

                    phase2_seen, regions_count, failed_count, trace_snapshot = _run_phase2_window(
                        phase2_inputs, options, picklable_options_p2, user_workers,
                        json_writer, all_modules, options_layer, data_writer, gbk_handle,
                        timings_by_record, lightweight_records, record_summaries,
                        detection_module_names, phase2_seen, regions_count, failed_count,
                        trace_snapshot, window_index, phase2_window_size,
                    )
                    gc.collect()

            if window_index == 0:
                logging.info("No records with regions found, skipping Phase 2")

            json_writer.finalize(timings_by_record)
        finally:
            if data_writer is not None:
                data_writer.finalize()
            elif html_enabled:
                data_file_path = os.path.join(options.output_dir, "regions_data.js")
                with open(data_file_path, "w", encoding="utf-8") as empty_data_handle:
                    StreamingRegionDataWriter(empty_data_handle).finalize()
            if data_handle is not None:
                data_handle.close()
            if gbk_handle is not None:
                gbk_handle.close()

    # Report failures
    if failed_count:
        logging.warning("Streaming pipeline: %d of %d records failed",
                        failed_count, accepted_count)
    if total_processed == 0:
        logging.error("All %d records failed, aborting", accepted_count)
        return 1

    skipped_count = total_without_regions if skip_without_regions else 0

    logging.info("Streaming complete: %d with regions, %d without, %d failed",
                 regions_count, total_without_regions, failed_count)

    # HTML finalization (overview page + dashboard) — uses lightweight data only
    html_start = time.time()
    if html_enabled:
        taxonomy_mapping: Dict[str, str] = {}
        if options.html_taxonomy:
            taxonomy_mapping = parse_taxonomy_file(options.html_taxonomy)

        finalize_streaming_html_output(
            lightweight_records, record_summaries,
            options, all_modules, taxonomy_mapping=taxonomy_mapping,
            skipped_record_count=skipped_count,
        )
        trace_snapshot = _log_parent_memory(
            "html-finalize-complete",
            options,
            extra={
                "regions_count": regions_count,
                "stored_record_summaries": len(record_summaries),
                "stored_lightweight_records": len(lightweight_records),
                "skipped": skipped_count,
            },
            trace_snapshot=trace_snapshot,
        )
        html_duration = (time.time() - html_start) / max(regions_count, 1)
        for val in timings_by_record.values():
            val[html.__name__] = html_duration

    # GenBank output for records WITHOUT regions was handled during Phase 2;
    # only zip remains
    base_filename = canonical_base_filename(input_basename, options.output_dir, options)
    zipfile = base_filename + ".zip"
    if os.path.exists(zipfile):
        os.remove(zipfile)
    if options.zip_output:
        logging.debug("Zipping output to '%s'", zipfile)
        with tempfile.NamedTemporaryFile(prefix="as_zip_tmp", suffix=".zip") as temp:
            shutil.make_archive(temp.name.replace(".zip", ""), "zip", root_dir=options.output_dir)
            shutil.copy(temp.name, zipfile)
            os.chmod(zipfile, 0o644)
        assert os.path.exists(zipfile)

    running_time = datetime.now() - start_time
    trace_snapshot = _log_parent_memory(
        "streaming-complete",
        options,
        extra={
            "processed": total_processed,
            "regions_count": regions_count,
            "without_regions": total_without_regions,
            "failed": failed_count,
        },
        trace_snapshot=trace_snapshot,
    )

    if options.debug:
        log_module_runtimes(timings_by_record)

    logging.debug("antiSMASH calculation finished at %s; runtime: %s",
                  datetime.now().strftime("%Y-%m-%d %H:%M:%S"), str(running_time))

    logging.info("antiSMASH status: SUCCESS")
    return 0


def _create_picklable_options(options: ConfigType) -> ConfigType:
    """Create a copy of the options object without unpicklable module objects.

    The all_enabled_modules attribute contains Python module objects which
    cannot be pickled for multiprocessing. This creates a copy without them.

    Arguments:
        options: the original options object

    Returns:
        a copy of the options object without module objects
    """
    picklable_options = copy.copy(options)

    # Remove the all_enabled_modules attribute, which contains module objects
    if hasattr(picklable_options, 'all_enabled_modules'):
        delattr(picklable_options, 'all_enabled_modules')

    return picklable_options


METASMASH_BUILD = "ms-build-7.5"


def _run_antismash(sequence_file: Optional[str], options: ConfigType) -> int:
    """ The real run_antismash, assumes logging is set up around it """
    logging.info("antiSMASH version: %s", options.version)
    logging.info("MetaSMASH build: %s", METASMASH_BUILD)
    _log_found_executables(options)

    if options.list_plugins:
        list_plugins()
        return 0

    modules = get_all_modules()
    options.all_enabled_modules = _get_all_enabled_modules(modules, options)

    if options.check_prereqs_only:
        try:
            check_prerequisites(modules, options)
        except RuntimeError:
            print("Some module prerequisites not satisfied")
            return 1
        print("All prerequisites satisfied")
        return 0

    check_prerequisites(options.all_enabled_modules, options)

    # start up profiling if relevant
    if options.profile:
        profiler = cProfile.Profile()
        profiler.enable()

    # ensure the provided options are valid
    if not verify_options(options, options.all_enabled_modules):
        return 1

    # check that at least one module will run
    if not options.all_enabled_modules:
        raise ValueError("No detection or analysis modules enabled")

    # resolve --workers: default to --cpus for backward compatibility
    if options.workers < 1:
        update_config({"workers": options.cpus})
    else:
        update_config({"workers": min(options.workers, options.cpus)})

    # determine if streaming mode should be used
    use_streaming, prefetched_metadata = should_use_streaming(options, sequence_file)

    if use_streaming:
        logging.info("Using streaming record processing pipeline")
        result = _run_antismash_streaming(sequence_file, options, prefetched_metadata)
        # save profiling data
        if options.profile:
            profiler.disable()
            write_profiling_results(profiler, os.path.join(options.output_dir,
                                                           "profiling_results"))
        return result

    start_time = datetime.now()

    results = read_data(sequence_file, options)

    # reset module timings
    results.timings_by_record.clear()

    prepare_output_directory(options.output_dir, sequence_file or options.reuse_results)

    results.records = record_processing.pre_process_sequences(results.records, options,
                                                              cast(AntismashModule, genefinding))

    # Create picklable options (removes all_enabled_modules which contains unpicklable module objects)
    from antismash.common.subprocessing import parallel_function

    picklable_options = _create_picklable_options(options)

    # Build list of non-skipped records for parallel processing
    record_and_results_list = [(record, mod_results)
                               for record, mod_results in zip(results.records, results.results)
                               if not record.skip]

    if record_and_results_list:
        _preload_pfam_caches(options)
        _preload_analysis_caches(options)
        gc.freeze()

        # Create a mapping of record IDs to their indices
        record_id_to_index = {record.id: i for i, record in enumerate(results.records)}

        # Phase 1: Parallel detection
        logging.info("Starting parallel processing of %d records for detection", len(record_and_results_list))
        detection_output = parallel_function(
            process_record_detection,
            [(rec_res, picklable_options) for rec_res in record_and_results_list],
            cpus=options.workers,
        )

        # Update records/results from detection output
        for record, mod_results, timings in detection_output:
            j = record_id_to_index[record.id]
            results.records[j] = record
            results.results[j] = mod_results
            if timings:
                results.timings_by_record[record.id] = timings

        # Phase 2: Parallel analysis (only records with regions)
        analysis_list = [(record, results.results[record_id_to_index[record.id]])
                         for record, _, _ in detection_output
                         if not record.skip and record.get_regions()]

        if analysis_list:
            logging.info("Starting parallel processing of %d records for analysis", len(analysis_list))
            analysis_output = parallel_function(
                process_record_analysis,
                [(rec_res, picklable_options) for rec_res in analysis_list],
                cpus=options.workers,
            )

            # Update records/results from analysis output
            for record, mod_results, timings in analysis_output:
                j = record_id_to_index[record.id]
                results.records[j] = record
                results.results[j] = mod_results
                if timings:
                    if record.id in results.timings_by_record:
                        results.timings_by_record[record.id].update(timings)
                    else:
                        results.timings_by_record[record.id] = timings

    # Write results (apply --skip-records-without-regions for JSON)
    logging.info("Writing results")
    filtered_records, filtered_results_list, _skip_count, _skip_ids = \
        _filter_records_for_output(results.records, results.results, options)
    json_results = serialiser.AntismashResults(
        results.input_file, filtered_records, filtered_results_list,
        results.version, timings=results.timings_by_record, taxon=results.taxon)
    json_filename = canonical_base_filename(results.input_file, options.output_dir, options)
    json_filename += ".json"
    logging.debug("Writing json results to '%s'", json_filename)
    json_results.write_to_file(json_filename)

    # now that the json is out of the way, annotate the record
    # otherwise we could double annotate some areas
    annotate_records(results)

    # create relevant output files (uses original results for HTML, filtered for GBK)
    write_outputs(results, options)

    # save profiling data
    if options.profile:
        profiler.disable()
        write_profiling_results(profiler, os.path.join(options.output_dir,
                                                       "profiling_results"))

    running_time = datetime.now() - start_time

    # display module runtimes before total time
    if options.debug:
        log_module_runtimes(results.timings_by_record)

    logging.debug("antiSMASH calculation finished at %s; runtime: %s",
                  datetime.now().strftime("%Y-%m-%d %H:%M:%S"), str(running_time))

    logging.info("antiSMASH status: SUCCESS")
    return 0


def _log_found_executables(options: ConfigType) -> None:
    for binary, path in vars(options.executables).items():
        version = ""
        version_getter = getattr(subprocessing, f"run_{binary}_version", None)
        if callable(version_getter):
            version = f" ({version_getter()})"
        logging.info("%s using executable: %s%s", binary, path, version)
