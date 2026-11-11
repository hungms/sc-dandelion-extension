# dandelion_extension

Scripts for running [sc-dandelion](https://github.com/zktuong/dandelion) BCR/TCR
preprocessing and contig QC inside the official dandelion Singularity/Apptainer container.

| Script | What it does |
|---|---|
| `dandelion_preprocess.py` | A modified version of the container's own `/share/dandelion_preprocess.py` (format contigs, re-annotate with IgBLAST, TIgGER allele reassignment, isotype assignment, mutation counts). See [Changes from upstream](#changes-from-upstream). |
| `dandelion_qc.py` | A new step. It merges the per-sample results, filters contigs, quantifies mutations per locus, calls clonotypes and writes AIRR / h5ddl / metadata outputs. |

This repo holds **code only**. It has no container images and no cluster-specific settings.

## How the setup works

```
┌──────────────────────────┐   bind-mount (-B)   ┌───────────────────────────────┐
│ this repo (git clone)    │ ──────────────────► │ dandelion container (.sif)    │
│  dandelion_preprocess.py │                     │  python, dandelion, IgBLAST,  │
│  dandelion_qc.py         │                     │  reference DBs (unchanged)    │
└──────────────────────────┘                     └───────────────────────────────┘
             ▲                                                  ▲
             │ path to clone                                    │ path to .sif
┌────────────┴──────────────────────────────────────────────────┴──────────────┐
│ your site wrapper (NOT in this repo): scheduler script, module loads,        │
│ data locations, per-dataset options                                          │
└──────────────────────────────────────────────────────────────────────────────┘
```

1. **Container**: the stock dandelion image supplies all the software and reference
   databases. Download it once and keep it wherever you like. It is never modified or rebuilt.
2. **Code (this repo)**: clone it once. At run time the scripts are bind-mounted into the
   container and run with the container's Python. Pulling an update or editing a script
   therefore takes effect straight away: no container rebuild, and no copying scripts next to
   each dataset.
3. **Site wrapper (yours)**: a small job script for your scheduler that sets resources, loads
   Singularity/Apptainer, and calls the two scripts with your dataset's options. Keep it
   outside this repo (for example in a separate repo of pipelines for your cluster), so that
   cluster details never get mixed into the shared code. A minimal Slurm example is
   [below](#running-on-a-cluster).

## Requirements

- Singularity ≥ 3.x or Apptainer
- The sc-dandelion container. Tested with **sc-dandelion 0.5.4** (image pulled 2026-03-06):
  ```bash
  singularity pull library://kt16/default/sc-dandelion:latest   # -> sc-dandelion_latest.sif
  ```
  `latest` changes over time. Keep the `.sif` you used for a project so results stay reproducible.
- Cell Ranger VDJ outputs: `filtered_contig.fasta` + `filtered_contig_annotations.csv`
  (or the `all_` equivalents) for each sample.

## Install

```bash
git clone https://github.com/hungms/dandelion_extension.git
```

## Input layout

Both scripts run from a single working directory that holds one sub-folder per sample:

```
work_dir/
├── sample.csv                          # sample-level metadata (preprocess)
├── cell.csv                            # cell-level metadata (QC)
├── SAMPLE_1/
│   ├── filtered_contig.fasta
│   └── filtered_contig_annotations.csv
└── SAMPLE_2/
    └── ...
```

**`sample.csv`** (passed with `--meta` to preprocess): has a header row. The first column is
the sample ID, which must match the folder name. Optional columns:
- `prefix` / `suffix`: text added to barcodes
- `individual`: groups samples for TIgGER (otherwise all samples are analysed together)

**`cell.csv`** (passed with `--meta` to QC): has a header row, and the first column (index) is the
cell barcode, written the same way as in the preprocessed output (e.g. `SAMPLE_1_AAAC...`). It
**must have a `project` column** whose values are the sample folder names. These values tell
the script which `SAMPLE/dandelion/<prefix>_contig_dandelion.tsv` files to read.

## Usage

```bash
CODE=/path/to/dandelion_extension
SIF=/path/to/sc-dandelion_latest.sif
cd /path/to/work_dir

# 1. Preprocess every sample folder
singularity run -B "$PWD,$CODE" "$SIF" python "$CODE/dandelion_preprocess.py" \
    --org mouse --chain IG --db imgt --file_prefix filtered \
    --meta sample.csv --keep_trailing_hyphen_number

# 2. Merge, filter contigs and call clonotypes
singularity run -B "$PWD,$CODE" "$SIF" python "$CODE/dandelion_qc.py" \
    --meta cell.csv --chain IG --file_prefix filtered --size 4
```

`-B` must include the working directory and the repo clone. Add any other directory that your
inputs come from. To list all options, run
`singularity run "$SIF" python "$CODE/<script>.py" --help`.

**Troubleshooting:** `RuntimeError: cannot cache function ... no locator available` means
numba (imported by scanpy) cannot write its cache inside the read-only container. Point it at
a writable directory:
```bash
export SINGULARITYENV_NUMBA_CACHE_DIR=/tmp   # APPTAINERENV_NUMBA_CACHE_DIR for Apptainer
```

### Outputs

- Preprocess: `SAMPLE_*/dandelion/<prefix>_contig_dandelion.tsv` (plus TIgGER/isotype plots)
- QC, in `work_dir/contig_qc/`:
  - `dandelion_airr.tsv`: filtered contigs in AIRR format
  - `dandelion_data.h5ddl`: dandelion object
  - `dandelion_meta.csv`: per-cell metadata with dandelion annotations

### `dandelion_preprocess.py` options

| Option | Default | Meaning |
|---|---|---|
| `--meta` | none | Sample metadata CSV (see above). Without it, every non-hidden sub-folder is treated as a sample. |
| `--org` | `human` | `human` or `mouse` |
| `--chain` | `IG` | `IG` or `TR` (TR stops after re-annotation) |
| `--db` | `imgt` | `imgt` or `ogrdb` |
| `--strain` | all | Mouse strain (only with `--db ogrdb`) |
| `--file_prefix` | `all` | Uses `<prefix>_contig.fasta` / `<prefix>_contig_annotations.csv` |
| `--sep` | `_` | Separator between barcode and prefix/suffix |
| `--flavour` | `strict` | IgBLAST mode: `strict` (applies e-value/penalty cut-offs) or `original` |
| `--filter_to_high_confidence` | off | Keep only contigs with `high_confidence == True` |
| `--keep_trailing_hyphen_number` | off | Keep the `-1` suffix on barcodes |
| `--skip_format_header` | off | Skip contig header formatting |
| `--skip_tigger` | off | Skip TIgGER allele reassignment |
| `--skip_reassign_dj` | off | Skip blastn D/J reassignment (strict flavour) |
| `--skip_correct_c` | off | Has no effect in this version: C-call correction is always off (see below) |
| `--clean_output` | off | Delete intermediate files |

### `dandelion_qc.py` options

| Option | Default | Meaning |
|---|---|---|
| `--meta` | required | Cell metadata CSV with a `project` column (see above) |
| `--chain` | `IG` | `IG` or `TR`. BCR clones are called at 85% CDR3 amino-acid identity; TCR clones by exact CDR3 nucleotide sequence. |
| `--file_prefix` | `filtered` | Reads `SAMPLE/dandelion/<prefix>_contig_dandelion.tsv` |
| `--skip_correct_c` | off | Replace `c_call` with Cell Ranger's `c_call_10x` |
| `--skip_find_clones` | off | Do not call clonotypes |
| `--size` | `4` | `max_size` for the clone-size categories |

Contig filtering in QC keeps the highest-UMI contig (UMI fold-change cut-off 2), removes extra
VDJ chains and missing contigs, and keeps non-productive contigs. Mutation counts for IG are
split by locus.

## Running on a cluster

Wrap the commands above in a job script for your scheduler and keep that script outside this
repo. A minimal Slurm example:

```bash
#!/bin/bash
#SBATCH --job-name=dandelion_preprocess
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=2-00:00
set -euo pipefail
module load singularity        # or apptainer: use whatever your site provides

CODE=/path/to/dandelion_extension
SIF=/path/to/sc-dandelion_latest.sif
cd /path/to/work_dir

singularity run -B "$PWD,$CODE" "$SIF" python "$CODE/dandelion_preprocess.py" \
    --org mouse --chain IG --meta sample.csv --file_prefix filtered
```

Submit QC as a second job that waits for preprocessing to succeed:

```bash
jid=$(sbatch --parsable preprocess.sh)
sbatch --dependency=afterok:"$jid" qc.sh
```

To run a fixed version for a production analysis, check out a tag in the clone
(`git -C "$CODE" checkout v0.1.0`) and record it with your results.

## Changes from upstream

Compared with `/share/dandelion_preprocess.py` in the sc-dandelion 0.5.4 container:

- `filename_prefix` is passed to every dandelion function as a per-sample list instead of a
  single string.
- **C-call correction is disabled** (`correct_c_call=False` at `assign_isotypes`) whatever
  `--skip_correct_c` is set to. A mouse IGHG2A/B/C correction dictionary is defined but not used.
- Sample names come straight from the metadata index, without converting them to strings.
- This script was based on an earlier upstream version, so it does **not** include upstream's
  later safety check that stops `--clean_output` from deleting sample folders when an
  `individual` value has the same name as a sample. Avoid using `--clean_output` with
  overlapping names.

`dandelion_qc.py` is new and has no upstream counterpart.

## Licence and credit

Based on [sc-dandelion](https://github.com/zktuong/dandelion) by Zewen Kelvin Tuong and
contributors, which is licensed under the GNU Affero General Public License v3 or later. This
repo is distributed under the same licence (see `LICENSE`). If you use it, please cite
dandelion:

> Suo C, Polanski K, Dann E, et al. Dandelion uses the single-cell adaptive immune receptor
> repertoire to explore lymphocyte developmental origins. *Nature Biotechnology* 42, 40–51 (2024).
