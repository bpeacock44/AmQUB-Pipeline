#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path

from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# ----------------------------
# CONSTANTS
# ----------------------------

PH_TARGETS = {
    "H_brev_a1","H_ulicic1","H_tortuo1","H_lunata1",
    "H_juliae1","H_aff_j1","H_subfus1","H_inflat1",
    "H_fusisp1","H_multig1","H_brev_s1","H_fagi1",
    "H_fusari1","ARF1","DoUCR501","GloDo1",
    "HsImV251","HsImV271","StM1","H_roti1"
}

OUTGROUP = "Lecophag1"

# ----------------------------
# FASTA PROCESSING
# ----------------------------

def combine_fastas(ref_file, input_file, output_file):
    with open(output_file, "w") as out:
        with open(ref_file) as r:
            out.write(r.read())
        with open(input_file) as i:
            out.write(i.read()) 

def write_ph_pph_fastas(fasta_file, ph_set, pph_set, outdir, base):
    ph_out = outdir / f"{base}_PH.fa"
    pph_out = outdir / f"{base}_PPH.fa"

    ph_records = []
    pph_records = []

    for rec in SeqIO.parse(fasta_file, "fasta"):
        name = rec.id

        if name in ph_set:
            ph_records.append(rec)
        elif name in pph_set:
            pph_records.append(rec)

    SeqIO.write(ph_records, ph_out, "fasta")
    SeqIO.write(pph_records, pph_out, "fasta")

# ----------------------------
# ALIGNMENT
# ----------------------------

def run_mafft(input_fa, output_aln, threads=8):
    cmd = [
        "mafft",
        "--thread", str(threads),
        "--auto",
        input_fa
    ]
    with open(output_aln, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)

# ----------------------------
# TREE BUILDING
# ----------------------------

def build_nj_tree(alignment_file):
    aln = AlignIO.read(alignment_file, "fasta")
    calc = DistanceCalculator("identity")
    dm = calc.get_distance(aln)
    constructor = DistanceTreeConstructor()
    return constructor.nj(dm)

# ----------------------------
# PH CLASSIFICATION
# ----------------------------

def classify_ph(tree):
    leaves = list(tree.get_terminals())
    ph_leaves = [t for t in leaves if t.name in PH_TARGETS]
    if not ph_leaves:
        return None, None
    clade = tree.common_ancestor(ph_leaves)
    ph = {t.name for t in clade.get_terminals()}
    all_taxa = {t.name for t in leaves}
    pph = all_taxa - ph
    return ph, pph

# ----------------------------
# MAIN PIPELINE PER FILE
# ----------------------------

def process_file(fa_file, ref_file, outdir, threads):
    base = Path(fa_file).stem
    print(f"Processing {base}")

    combined = outdir / f"{base}.ready.fa"
    aln = outdir / f"{base}.aln"

    combine_fastas(ref_file, fa_file, combined)

    if not combined.exists():
        raise RuntimeError(f"Combined FASTA missing: {combined}")

    run_mafft(combined, aln, threads=threads)

    combined.unlink(missing_ok=True)

    # Tree
    tree = build_nj_tree(aln)

    try:
        tree.root_with_outgroup(OUTGROUP)
    except Exception:
        print(f"  Warning: {OUTGROUP} not found in {base}")

    nj_path = outdir / f"{base}.nj.nwk"
    Phylo.write(tree, nj_path, "newick")

    # PH classification
    ph, pph = classify_ph(tree)

    if ph is None:
        print(f"  No PH targets found")
        return

    ph_filtered = ph - PH_TARGETS
    pph_filtered = pph - PH_TARGETS - {OUTGROUP}

    with open(outdir / f"{base}_PH.txt", "w") as f:
        f.write("\n".join(sorted(ph_filtered)))

    with open(outdir / f"{base}_PPH.txt", "w") as f:
        f.write("\n".join(sorted(pph_filtered)))

    write_ph_pph_fastas(
        fa_file,
        ph_filtered,
        pph_filtered,
        outdir,
        base
    )

# ----------------------------
# ENTRYPOINT
# ----------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("base_dir", help="Directory containing .fa files")
    parser.add_argument("--ref", required=True, help="Reference FASTA")
    parser.add_argument("--threads", type=int, default=8)

    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    outdir = base_dir

    fasta_files = sorted(base_dir.glob("*.fa"))

    for fa in fasta_files:
        process_file(fa, args.ref, outdir, args.threads)

    print("Processing complete.")

if __name__ == "__main__":
    main()