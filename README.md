# MetaSMASH

A scalable fork of [antiSMASH](https://github.com/antismash/antismash) for metagenome-scale biosynthetic gene cluster (BGC) mining.

## What's different

MetaSMASH adds a streaming pipeline on top of antiSMASH that scales to millions of input records with bounded memory usage:

- **Streaming pipeline** (`--streaming auto|on|off`) — two-phase processing (detection, then analysis) that processes records one at a time instead of loading all into memory. Auto-enabled for inputs with >10 records.
- **Incremental output** — JSON and GenBank files are written per-record as they complete, so results are available even if a run is interrupted.
- **Lazy parallel execution** — worker results are yielded one at a time instead of collected, enabling constant-memory parallelism.
- **Cache preloading with fork CoW** — PFAM and analysis databases are loaded before forking so worker processes share read-only memory via copy-on-write.
- **Metagenome dashboard** — summary view across all records in the HTML output.
- **Record filtering** — `--output-skip-records-without-regions` (on by default) omits records without detected BGC regions from the output.

The classic (non-streaming) antiSMASH pipeline is fully preserved and can be used with `--streaming off`.

## Installation

Requires Python >= 3.11. Not on PyPI — install from source.

```bash
# Create a conda environment with required external tools
conda create -n metasmash -c bioconda -c conda-forge \
    "python>=3.11" blast diamond hmmer hmmer2 prodigal

conda activate metasmash

# Clone and install MetaSMASH
git clone https://github.com/canerbagci/metasmash.git
cd metasmash
pip install .

# Download databases (same as upstream antiSMASH)
download-antismash-databases
```

## Usage

```bash
# Basic usage (streaming auto-enabled for large inputs)
metasmash input.fasta

# Force streaming mode
metasmash --streaming on input.fasta

# Disable streaming (classic antiSMASH behaviour)
metasmash --streaming off input.fasta
```

### New options

| Option | Default | Description |
|--------|---------|-------------|
| `--streaming {auto,on,off}` | `auto` | `auto` enables streaming for >10 records, `on` forces it, `off` uses classic pipeline |
| `--workers W` | same as `--cpus` | Number of parallel worker processes. Each worker gets `cpus/workers` threads. Lower values reduce peak memory. |
| `--output-skip-records-without-regions` / `--no-output-skip-records-without-regions` | on | Omit records without detected BGC regions from JSON/GBK output |

All standard antiSMASH options (e.g. `--cpus`, `--genefinding-tool`, `--cb-knownclusters`, `--minimal`) work as usual. Run `metasmash --help` for the full list, or see the [upstream documentation](https://docs.antismash.secondarymetabolites.org/).

### Parallelism and memory tuning

`--cpus` and `--workers` together control the trade-off between throughput and memory usage:

- **`--cpus N`** (default: all cores) — total CPU budget for the run.
- **`--workers W`** (default: same as `--cpus`) — number of parallel worker processes. Each worker gets `cpus / workers` threads for internal tools (hmmsearch, diamond, blastp, etc.).

In streaming mode, the two phases use workers differently:

| Phase | Workers | Threads per worker | Purpose |
|-------|---------|-------------------|---------|
| Phase 1 (detection) | `cpus` | 1 | Many lightweight workers scanning records in parallel |
| Phase 2 (analysis) | `workers` | `cpus / workers` | Fewer workers with more threads for heavier analysis modules |

**Tuning guidance:**

- **Fewer workers = less memory.** Each worker loads one record into memory, so fewer workers means fewer records resident simultaneously. Each worker compensates with more threads for I/O-heavy tools.
- **More workers = higher throughput, higher peak memory.** More records processed in parallel, but each holds its own copy of per-record data structures.
- **Example:** `metasmash --cpus 32 --workers 8 input.fasta` runs 8 parallel records, each using 4 threads internally — a good balance for large metagenomes on a 32-core machine.
- **Default (workers = cpus):** maximises record-level parallelism with single-threaded tools. Fine for small-to-medium inputs; consider lowering `--workers` for very large metagenomes to keep memory bounded.

## Upstream antiSMASH

For general antiSMASH usage, module documentation, citation information, and the web server, see the upstream project:

- **Repository**: <https://github.com/antismash/antismash>
- **Documentation**: <https://docs.antismash.secondarymetabolites.org/>
- **Citations**: <http://antismash.secondarymetabolites.org/#!/about>

## License

MetaSMASH is available under the GNU Affero General Public License v3.0 or later, same as upstream antiSMASH. See [`LICENSE.txt`](LICENSE.txt) for details.
