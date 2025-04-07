from ete3 import Tree
import argparse
import os

def load_tree(tree_path, format_val, quoted_names):
    print(f"Loading tree from: {tree_path}")
    print(f"Using format={format_val}, quoted_node_names={quoted_names}")
    t = Tree(tree_path, format=format_val, quoted_node_names=quoted_names)
    print(f"Tree loaded: {len(t)} leaves, {len(t.get_descendants())} total nodes")
    return t

def prune_tree(tree, taxa_list):
    print(f"Pruning tree to {len(taxa_list)} taxa...")
    tree.prune(taxa_list, preserve_branch_length=True)
    print(f"Pruned tree has {len(tree)} leaves")
    return tree

def save_tree(tree, out_path):
    tree.write(format=1, outfile=out_path)
    print(f"Tree written to: {out_path}")

def main(args):
    tree = load_tree(args.tree, args.format, args.quoted_node_names)

    if args.prune_to:
        with open(args.prune_to, "r") as f:
            taxa_list = [line.strip() for line in f if line.strip()]
        tree = prune_tree(tree, taxa_list)

    save_tree(tree, args.out)

    if args.display:
        print("Launching interactive tree viewer...")
        tree.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="MK's script for Dharmik's phylogenetic tree bias -- newick format in ETE Toolkit.")
    parser.add_argument("--tree", required=True, help="Path -- Newick tree")
    parser.add_argument("--out", required=True, help="Output path for (pruned) Newick tree")
    parser.add_argument("--prune_to", help="Path to text file with taxa to retain (one per line)")
    parser.add_argument("--display", action="store_true", help="Display tree -- (GUI required)")
    parser.add_argument("--format", type=int, default=1, help="Tree parsing format code (default: 1)")
    parser.add_argument("--quoted_node_names", action="store_true", help="Use this if tree tip/node names are quoted")

    args = parser.parse_args()
    main(args)
