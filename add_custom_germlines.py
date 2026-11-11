#!/opt/conda/envs/sc-dandelion-container/bin/python
"""Add custom V/D/J genes to a copy of the sc-dandelion reference database.

Run inside the dandelion container on a *copy* of /share/database, then bind that copy over
/share/database when running the pipeline (see README, "Adding custom germline genes").

    python add_custom_germlines.py <db_dir> <new_genes.fasta> [--org mouse]

<new_genes.fasta> holds ungapped sequences. Name each gene in IMGT style so that Change-O and
TIgGER can parse it (e.g. IGHV0-HyHel10*01). V genes also need the closest germline allele,
whose IMGT gap pattern is copied onto the new sequence:

    >IGHV0-HyHel10*01 template=IGHV3-8*02
    GACGTGCAGCTTCAGGAG...
    >IGHD0-HyHel10*01
    AACTGGGACGGTGATTAC

For each gene this script:
  1. appends the IMGT-gapped sequence to germlines/imgt/<org>/vdj/imgt_<org>_<IGHV|IGHD|...>.fasta
     (read by MakeDb, create_germlines and TIgGER);
  2. appends the ungapped sequence to igblast/fasta/imgt_<org>_<ig|tr>_<v|d|j>.fasta;
  3. rebuilds the matching BLAST database in igblast/database/ (searched by IgBLAST and by
     dandelion's blastn D/J re-alignment).
Genes that are already present are skipped, so running it twice is safe.
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Change-O's allele pattern (changeo/Gene.py); names that fail it lose their gene call.
ALLELE_RE = re.compile(r"^(IG[HKL]|TR[ABGD])([VDJ])[A-R0-9]+[-/\w]*[-*][.\w]+$")
MIN_IDENTITY = 0.85


def read_fasta(path):
    """Return {header: sequence} for a FASTA file (multi-line sequences allowed)."""
    seqs, header = {}, None
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if line.startswith(">"):
            header = line[1:]
            seqs[header] = ""
        elif header is not None:
            seqs[header] += line
    return seqs


def gap_like(seq, template):
    """Copy the IMGT gaps ('.') of a gapped template onto an ungapped sequence."""
    out, i = [], 0
    for c in template:
        if i == len(seq):
            break
        if c == ".":
            out.append(".")
        else:
            out.append(seq[i])
            i += 1
    out.append(seq[i:])  # any 3' bases beyond the template
    ungapped = template.replace(".", "")
    n = min(len(seq), len(ungapped))
    identity = sum(a == b for a, b in zip(seq[:n], ungapped[:n])) / n
    return "".join(out), identity


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("db_dir", type=Path, help="copy of the container's /share/database")
    parser.add_argument("new_genes", type=Path, help="FASTA of ungapped custom genes")
    parser.add_argument("--org", default="mouse")
    args = parser.parse_args()

    vdj_dir = args.db_dir / "germlines" / "imgt" / args.org / "vdj"
    igblast_dir = args.db_dir / "igblast"
    makeblastdb = shutil.which("makeblastdb") or str(Path(sys.executable).parent / "makeblastdb")
    to_rebuild = set()

    for header, seq in read_fasta(args.new_genes).items():
        name, *opts = header.split()
        opts = dict(o.split("=", 1) for o in opts)
        seq = seq.upper().replace(".", "")
        m = ALLELE_RE.match(name)
        if not m:
            sys.exit(f"{name}: not an IMGT-style allele name (e.g. IGHV0-HyHel10*01).")
        locus, segment = m.group(1), m.group(2)

        # 1. gapped germline (source of truth)
        germ_file = vdj_dir / f"imgt_{args.org}_{locus}{segment}.fasta"
        germ = read_fasta(germ_file)
        if any(h.split("|")[1] == name for h in germ):
            print(f"{name}: already in {germ_file.name}, skipped")
        else:
            gapped = seq
            if segment == "V":
                if "template" not in opts:
                    sys.exit(f"{name}: V genes need 'template=<closest allele>' in the header.")
                tmpl = [s for h, s in germ.items() if h.split("|")[1] == opts["template"]]
                if not tmpl:
                    sys.exit(f"{name}: template {opts['template']} not found in {germ_file.name}.")
                gapped, identity = gap_like(seq, tmpl[0].upper())
                if identity < MIN_IDENTITY:
                    sys.exit(
                        f"{name}: only {identity:.0%} identical to {opts['template']} position by "
                        "position (indels?). Gap it with IMGT/V-QUEST instead."
                    )
                print(f"{name}: gapped from {opts['template']} ({identity:.1%} identical)")
            n = len(seq)
            fields = ["custom", name, f"{args.org}_custom", "F", f"{segment}-REGION",
                      f"1..{n}", f"{n} nt", "1", " ", " ", " ", " ",
                      f"{n}+{len(gapped) - n}={len(gapped)}", " ", " ", ""]
            with open(germ_file, "a") as fh:
                fh.write(f"\n>{'|'.join(fields)}\n{gapped.lower()}\n")
            print(f"{name}: added to {germ_file.relative_to(args.db_dir)}")

        # 2. ungapped IgBLAST FASTA
        chain = "ig" if locus.startswith("IG") else "tr"
        ig_file = igblast_dir / "fasta" / f"imgt_{args.org}_{chain}_{segment.lower()}.fasta"
        if name in read_fasta(ig_file):
            print(f"{name}: already in {ig_file.name}, skipped")
        else:
            with open(ig_file, "a") as fh:
                fh.write(f"\n>{name}\n{seq}\n")
            print(f"{name}: added to {ig_file.relative_to(args.db_dir)}")
        to_rebuild.add(ig_file)

    # 3. rebuild the BLAST databases IgBLAST actually searches
    for ig_file in sorted(to_rebuild):
        # drop blank lines left by appending, which makeblastdb rejects
        text = "\n".join(l for l in ig_file.read_text().splitlines() if l.strip()) + "\n"
        ig_file.write_text(text)
        out = igblast_dir / "database" / ig_file.stem
        subprocess.run(
            [makeblastdb, "-parse_seqids", "-dbtype", "nucl", "-input_type", "fasta",
             "-in", str(ig_file), "-out", str(out)],
            check=True, stdout=subprocess.DEVNULL,
        )
        print(f"rebuilt BLAST database {out.relative_to(args.db_dir)}")


if __name__ == "__main__":
    main()
