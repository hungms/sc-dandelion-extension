#!/opt/conda/envs/sc-dandelion-container/bin/python
import argparse
import os
import shutil
import warnings

import dandelion as ddl
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from pathlib import Path
from scanpy import logging as logg

sc.settings.verbosity = 3

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--meta",
        help=(
            "Require cell-level metadata CSV file, header required"
            + 'Must include a column named "sample"'
            + 'Make sure each cell barcode has a prefix of "sample"_'
        ),
    )
    parser.add_argument(
        "--chain",
        type=str,
        default="IG",
        help=(
            'Indicate whether the data is TR or IG, as the preprocessing pipelines differ. Defaults to "IG".'
        ),
    )
    parser.add_argument(
        "--file_prefix",
        type=str,
        default="filtered",
        help=(
            "Which set of contig files to take for the folder. For a given "
            + 'PREFIX, will use PREFIX_contig_dandelion.tsv'
        ),
    )
    parser.add_argument(
        "--skip_correct_c",
        action="store_true",
        help=(
            'Whether to keep c_call from 10x'),
    )
    parser.add_argument(
        "--skip_find_clones",
        action="store_true",
        help=(
            'Whether to find clonotypes'),
    )
    parser.add_argument(
        "--size",
        type=int,
        default=4,
        help=('Max size per clonotype'),
    )
    args = parser.parse_args()
    print(args)

    # convert loci to lower case for compatibility, and ensure it's in TR/IG
    args.chain = args.chain.lower()
    if args.chain not in ["tr", "ig"]:
        raise ValueError("Chain must be TR or IG")
    return args


def main():
    """Filter Contigs & Clonotype Process."""
    logg.info("Software versions:\n")
    ddl.logging.print_header()
    # sponge up command line arguments to begin with
    args = parse_args()
    start = logg.info("\nBegin preprocessing\n")

    logg.info(
        "command line parameters:\n",
        deep=(
            f"\n"
            f"--------------------------------------------------------------\n"
            f"    --meta = {args.meta}\n"
            f"    --chain = {args.chain}\n"
            f"    --file_prefix = {args.file_prefix}\n"
            f"    --skip_correct_c = {args.skip_correct_c}\n"
            f"    --skip_find_clones = {args.skip_find_clones}\n"
            f"    --size = {args.size}\n"
            f"    beware TCR clones will be called by exact CDR3 nucleotide sequence, whilst BCR clones will called by 85% CDR3 amino acid similarity.\n"
            f": --------------------------------------------------------------\n"
        ),
    )

    # set TR/IG hamming distance threshold
    if args.chain == "tr":
        args.chain = "tr-ab"
        hamming = 1
        sequence = "junction"
    else:
        hamming = 0.85
        sequence = "junction_aa"

    # Step1 : read metadata & make anndata object
    meta = pd.read_csv(args.meta, index_col=0)
    counts = pd.DataFrame(np.zeros((meta.shape[0], 1)), columns=['x'])  # Corrected DataFrame
    counts.index.names = meta.index.names
    adata = sc.AnnData(X=counts, obs=meta[["project"]])

    # Step2 : set dandelion tsv files to read
    adatas = []
    samples = []
    column = adata.obs["project"].unique()
    for i in column:
        file = ''.join([i, '/dandelion/', args.file_prefix, '_contig_dandelion.tsv'])
        samples.append(file)
        mask = adata.obs["project"] == i
        adatas.append(adata[mask].copy())

    # Step3 : read tsv files transfer to adata
    if len(samples) != 0:
        vdjs = []
        for i in range(len(samples)):
            vdj = ddl.read_10x_airr(samples[i])
            # Step3 : skip_correct_c?
            if args.skip_correct_c:
                vdj.data["c_call"] = vdj.data["c_call_10x"]
                vdj.update_metadata()
            ddl.tl.transfer(adatas[i], vdj)
            vdjs.append(vdj)

        adatas = ad.concat(adatas)
        vdjs = ddl.concat(vdjs)

    # Step4 : remove poor quality contigs
    vdjs, adatas = ddl.pp.filter_contigs(
        vdjs,
        adatas,
        library_type=args.chain,
        filter_rna=False,
        filter_poorqualitycontig=False,
        keep_highest_umi=True,
        umi_foldchange_cutoff=2,
        filter_extra_vdj_chains=True,
        filter_extra_vj_chains=False,
        filter_missing=True,
        productive_only=False,
    )

    # Step5 : quantify mutations per locus
    if args.chain == "ig":
        ddl.pp.quantify_mutations(vdjs,split_locus=True,)
        ddl.pp.quantify_mutations(vdjs,frequency=True, split_locus=True,)

    # Step6 : find clones
    if not args.skip_find_clones:
        ddl.tl.find_clones(vdjs, identity=hamming, key=sequence)
        ddl.tl.clone_size(vdjs, max_size=args.size)

    # Step7 : save in airr format & save metadata
    vdjs.update_plus()
    vdjs.update_metadata()
    ddl.tl.transfer(adatas, vdjs)

    #contig_pass_meta = adatas.obs[adatas.obs['filter_contig'] == False]
    #contig_pass_meta.to_csv("contig_qc/dandelion_meta.csv")
    os.makedirs("contig_qc", exist_ok=True)
    adatas.obs.to_csv("contig_qc/dandelion_meta.csv")
    vdjs.write_airr("contig_qc/dandelion_airr.tsv")
    vdjs.write_h5ddl("contig_qc/dandelion_data.h5ddl")

    logg.info("Pre-processing finished.\n", time=start)


if __name__ == "__main__":
    main()
