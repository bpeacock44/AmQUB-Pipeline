#!/usr/bin/env python3
import argparse
import sys
import matplotlib.pyplot as plt
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


def find_tip(tree, query):
    terms = tree.get_terminals()
    for t in terms:
        if t.name == query:
            return t
    q = query.lower()
    matches = [t for t in terms if q in (t.name or "").lower()]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        print(f"Multiple matches for '{query}': {[m.name for m in matches]}",
              file=sys.stderr)
    return None


def aln_to_tree(aln_file, nwk_out, pdf_out, root_name="Lecophag1"):
    if not pdf_out:
        sys.exit("Error: pdf output path is empty")
    if not nwk_out:
        sys.exit("Error: newick output path is empty")

    aln = AlignIO.read(aln_file, "fasta")
    calc = DistanceCalculator("identity")
    dm = calc.get_distance(aln)
    tree = DistanceTreeConstructor().nj(dm)

    # Strip auto-generated "Inner###" labels on internal nodes
    # Strip auto-generated "Inner###" labels on internal nodes
    for clade in tree.get_nonterminals():
        clade.name = None

    target = find_tip(tree, root_name)
    if target is not None:
        tree.root_with_outgroup(target)
        print(f"Rooted at: {target.name}")
    else:
        print(f"Warning: could not find '{root_name}'. Available tips:",
              file=sys.stderr)
        for t in tree.get_terminals():
            print(f"  {t.name}", file=sys.stderr)

    Phylo.write(tree, nwk_out, "newick")

    def is_black(name):
        if not name:
            return False
        return name.startswith(("PH-m", "PH-O"))

    # Make blue labels uppercase
    for tip in tree.get_terminals():
        if not is_black(tip.name):
            tip.name = tip.name.upper()

    def label_color(name):
        if not name:
            return "black"
        return "black" if is_black(name) else "blue"

    n_leaves = tree.count_terminals()

    # Scale both dimensions to leaf count so big trees stay readable
    width = max(20, n_leaves * 0.08)
    height = max(10, n_leaves * 0.22)

    fig, ax = plt.subplots(figsize=(width, height))

    Phylo.draw(
        tree,
        axes=ax,
        label_colors=label_color,
        do_show=False
    )

    ax.set_title(f"NJ tree (rooted at {root_name})")

    plt.savefig(pdf_out, format="pdf", bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment")
    parser.add_argument("newick_output")
    parser.add_argument("pdf_output")
    parser.add_argument("--root", default="Lecophag1")
    args = parser.parse_args()
    aln_to_tree(args.alignment, args.newick_output, args.pdf_output, args.root)


if __name__ == "__main__":
    main()