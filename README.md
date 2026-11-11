# sc-dandelion-extension

Scripts for running [sc-dandelion](https://github.com/zktuong/dandelion) BCR/TCR
preprocessing and contig QC inside the official dandelion Singularity/Apptainer container.

| Script | What it does |
|---|---|
| `dandelion_preprocess.py` | A modified version of the container's own `/share/dandelion_preprocess.py` (format contigs, re-annotate with IgBLAST, TIgGER allele reassignment, isotype assignment, mutation counts). See [Changes from upstream](#changes-from-upstream). |
| `dandelion_preprocess_v0.3.2.py` | The same script adapted to the dandelion 0.3.2 API (see [Tested containers](#tested-containers)). |
| `dandelion_qc.py` | A new step. It merges the per-sample results, filters contigs, quantifies mutations per locus, calls clonotypes and writes AIRR / h5ddl / metadata outputs. |
| `add_custom_germlines.py` | Adds custom (e.g. transgenic) V/D/J genes to a copy of the reference database. See [Adding custom germline genes](#adding-custom-germline-genes-eg-swhel-hyhel10). |

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
- The sc-dandelion container. Download it, then rename the file after the dandelion version it
  contains:
  ```bash
  singularity pull library://kt16/default/sc-dandelion:latest   # -> sc-dandelion_latest.sif
  singularity exec sc-dandelion_latest.sif python -c "import dandelion; print(dandelion.__version__)"
  mv sc-dandelion_latest.sif sc-dandelion_latest_v<version>.sif
  ```
  `latest` changes over time. Keep the `.sif` you used for a project so results stay reproducible.
- Cell Ranger VDJ outputs: `filtered_contig.fasta` + `filtered_contig_annotations.csv`
  (or the `all_` equivalents) for each sample.
- See [Tested containers](#tested-containers) for which preprocess script goes with which
  container version.

## Install

```bash
git clone https://github.com/hungms/sc-dandelion-extension.git
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

### `sample.csv`

Passed with `--meta` to preprocess. One row per sample, with a header row.

| sample | prefix | individual |
|---|---|---|
| SAMPLE_1 | SAMPLE_1 | mouse_1 |
| SAMPLE_2 | SAMPLE_2 | mouse_1 |
| SAMPLE_3 | SAMPLE_3 | mouse_2 |

```csv
sample,prefix,individual
SAMPLE_1,SAMPLE_1,mouse_1
SAMPLE_2,SAMPLE_2,mouse_1
SAMPLE_3,SAMPLE_3,mouse_2
```

- **First column (required):** the sample ID. It must match the sample folder name.
- `prefix` or `suffix` (optional): text added to each barcode, joined with `--sep`
  (default `_`). Use one or the other, not both.
- `individual` (optional): groups samples from the same animal or donor for TIgGER. Without
  this column, all samples are genotyped together.

### `cell.csv`

Passed with `--meta` to QC. One row per cell, with a header row.

| *(index: cell barcode)* | project |
|---|---|
| SAMPLE_1_AAACCTGAGATGTTAG-1 | SAMPLE_1 |
| SAMPLE_1_AAACCTGAGTCCGGTC-1 | SAMPLE_1 |
| SAMPLE_2_AAACCTGCAAGCCGCT-1 | SAMPLE_2 |

```csv
"barcode","project"
"SAMPLE_1_AAACCTGAGATGTTAG-1","SAMPLE_1"
"SAMPLE_1_AAACCTGAGTCCGGTC-1","SAMPLE_1"
"SAMPLE_2_AAACCTGCAAGCCGCT-1","SAMPLE_2"
```

- **First column (index):** the cell barcode, written exactly as it appears in the
  preprocessed output: `<prefix><sep><barcode>`. The trailing `-1` is kept only if preprocess
  ran with `--keep_trailing_hyphen_number`. The header for this column can be blank (as in an
  `obs.to_csv()` export from AnnData).
- **`project` (required):** the sample folder name for each cell. It tells QC which
  `SAMPLE/dandelion/<prefix>_contig_dandelion.tsv` files to read.
- Extra columns are allowed. They're ignored.

## Usage

```bash
CODE=/path/to/sc-dandelion-extension
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

CODE=/path/to/sc-dandelion-extension
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

`dandelion_preprocess.py` is based on `/share/dandelion_preprocess.py` from the
**sc-dandelion 0.3.7** container (image pulled 2024-06-18). The only functional change is:

- **C-call correction is disabled** (`correct_c_call=False` at `assign_isotypes`) whatever
  `--skip_correct_c` is set to, so the mouse IGHG2A/B/C correction dictionary is not used.

The script runs unchanged in the 0.3.7 and 0.5.4 containers, because every dandelion function
it calls takes the same arguments in both versions. The 0.5.4 container ships a newer
`/share/dandelion_preprocess.py` that this repo has **not** picked up. That newer script
converts sample names to strings and adds a safety check so that `--clean_output` cannot delete
sample folders when an `individual` value has the same name as a sample. Avoid
`--clean_output` with overlapping names.

`dandelion_qc.py` is new and has no upstream counterpart.

## Licence and credit

Based on [sc-dandelion](https://github.com/zktuong/dandelion) by Zewen Kelvin Tuong and
contributors, which is licensed under the GNU Affero General Public License v3 or later. This
repo is distributed under the same licence (see `LICENSE`). If you use it, please cite
dandelion:

> Suo C, Polanski K, Dann E, et al. Dandelion uses the single-cell adaptive immune receptor
> repertoire to explore lymphocyte developmental origins. *Nature Biotechnology* 42, 40–51 (2024).

## Adding custom germline genes (e.g. SWHEL HyHEL10)

Transgenic or knock-in mice carry V(D)J genes that are not in IMGT. If those genes are missing,
the contigs are called against the closest endogenous gene, and every difference from it looks
like somatic hypermutation. `add_custom_germlines.py` adds the genes to a **copy** of the
container's reference database. You then bind that copy over `/share/database` at run time.
The container image itself is never changed.

The container keeps each gene in three places, and all three have to be updated together:

| Path inside `/share/database` | Contents | Read by |
|---|---|---|
| `germlines/imgt/<org>/vdj/imgt_<org>_IGHV.fasta` etc. | IMGT-gapped sequences (the source) | MakeDb, `create_germlines`, TIgGER, mutation counting |
| `igblast/fasta/imgt_<org>_ig_v.fasta` etc. | Ungapped copies | `makeblastdb` |
| `igblast/database/imgt_<org>_ig_v.n*` etc. | BLAST databases | IgBLAST, blastn D/J re-alignment |

Editing only `igblast/fasta/` has no effect, because IgBLAST searches the BLAST databases built
from it. `add_custom_germlines.py` updates all three.

**1. Write the new genes as ungapped FASTA.** For SWHEL, this is the HyHEL10 heavy-chain knock-in
(`swhel_hyhel10.fasta`):

```
>IGHV0-HyHel10*01 template=IGHV3-8*02
GACGTGCAGCTTCAGGAGTCAGGACCTAGCCTCGTGAAACCTTCTCAGACTCTGTCCCTCACCTGTTCTGTCACTGGCGACTCCATCACCAGTGATTATTGGAGCTGGATCCGGAAATTCCCAGGGAATAGACTTGAGTACATGGGGTACGTAAGCTACAGTGGTAGCACGTACTACAATCCATCTCTCAAAAGTCGAATCTCCATCACTCGAGACACATCCAAGAATCAGTACTATCTGGACCTGAATTCTGTGACTACTGAGGACACAGCCACATATTACTGTGCA
>IGHD0-HyHel10*01
AACTGGGACGGTGATTAC
```

- **Names** must follow the IMGT pattern `IG[HKL]` + `V/D/J` + family number + `-name*allele`,
  because Change-O and TIgGER read gene calls with that pattern. For example `IGHV0-HyHel10*01`
  works, but `HyHel10*01` or `IGHV-HyHel10*01` would lose their gene call. Family `0` marks the
  gene as custom.
- **V genes need `template=`**: the closest germline allele. The script copies that allele's
  IMGT gap pattern, which works as long as the new V has no insertions or deletions relative to
  it. It stops if identity is below 85%; in that case, gap the sequence with
  [IMGT/V-QUEST](https://www.imgt.org/IMGT_vquest/) instead.
- **Only add what isn't germline.** In SWHEL, the J region is an exact match to germline
  IGHJ3\*01 and the kappa chain is germline, so neither is added. The V ends at the codon after
  the second cysteine (`…TATTACTGTGCA`). The D covers the rest of the knock-in CDR3
  (IGHD4-1-like sequence plus N-additions), so the whole junction counts as germline.

**2. Build the custom database** (about 2 minutes):

```bash
SIF=/path/to/sc-dandelion_latest_v0.5.4.sif
CODE=/path/to/sc-dandelion-extension
DB=$PWD/database_v0.5.4_swhel

singularity exec -B "$PWD" "$SIF" cp -r /share/database "$DB"
singularity exec -B "$PWD,$CODE" "$SIF" python "$CODE/add_custom_germlines.py" \
    "$DB" swhel_hyhel10.fasta --org mouse
```

Running it again skips genes that are already present.

**3. Use it** by adding `$DB:/share/database` to `-B` on every `singularity run`:

```bash
singularity run -B "$PWD,$CODE,$DB:/share/database" "$SIF" python "$CODE/dandelion_preprocess.py" \
    --org mouse --chain IG --db imgt --file_prefix filtered --meta sample.csv
```

`dandelion_qc.py` reads the germlines stored in the preprocess output, so it doesn't need the
custom database. Bind it anyway to keep the two steps consistent.

## Tested containers

| Container file | dandelion | Python | Image pulled | Preprocess script |
|---|---|---|---|---|
| `sc-dandelion_latest_v0.3.2.sif` | 0.3.2 | 3.9.16 | 2023-11-14 | `dandelion_preprocess_v0.3.2.py` |
| `sc-dandelion_latest_v0.3.7.sif` | 0.3.7 | 3.11.9 | 2024-06-18 | `dandelion_preprocess.py` |
| `sc-dandelion_latest_v0.5.4.sif` | 0.5.4 | 3.12.9 | 2026-03-06 | `dandelion_preprocess.py` |

`dandelion_qc.py` is the same for all versions. It uses `ddl.pp.filter_contigs`, which is
deprecated in 0.5.x and scheduled for removal in 0.6.0.

`dandelion_preprocess_v0.3.2.py` is needed for 0.3.2 because that version has no OGRDB support
(so it takes only `--db imgt` and no `--strain`), and its `create_germlines` has a different
signature.

For SWHEL mice, use the stock 0.5.4 image with a custom database (see
[Adding custom germline genes](#adding-custom-germline-genes-eg-swhel-hyhel10)). An earlier
`sc-dandelion_latest_v0.5.4_swhel` setup edited only `igblast/fasta/` inside an unpacked
container and had no effect. Its `.sif` is identical to the stock image.
