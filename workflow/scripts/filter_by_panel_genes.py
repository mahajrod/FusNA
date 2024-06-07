#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="File with fusions called by Arriba. Default: stdin")
parser.add_argument("--g1", "--gene1_column_name", action="store", dest="gene1_column_name", default="#gene1",
                    help="Gene 1 column name. Default: #gene1")
parser.add_argument("--g2", "--gene2_column_name", action="store", dest="gene2_column_name", default="gene2",
                    help="Gene 2 column name. Default: gene2")
parser.add_argument("-c", "--control_genes", action="store", dest="control_genes", default=None,
                    help="File with IDs of control genes (one ID per row) used in the panel. Default: not set")
parser.add_argument("-t", "--target_genes", action="store", dest="target_genes", default=None,
                    help="File with IDs of target genes (One ID per row) used in panel. Default: not set")
parser.add_argument("--offtarget_output", action="store", dest="offtarget_output", required=True,
                    help="File with fusions involving offtarget genes.")
parser.add_argument("--control_output", action="store", dest="control_output", required=True,
                    help="File with fusions involving control genes.")
parser.add_argument("--target_output", action="store", dest="target_output", required=True,
                    help="File with fusions involving target genes.")
parser.add_argument("--all_output", action="store", dest="all_output", required=True,
                    help="File with all fusions labeled by gene type (target, control, offtarget).")
args = parser.parse_args()

panel_gene_dict = {}

for gene_type, gene_file in zip(["control", "target"],
                                [args.control_genes, args.target_genes]):
    if gene_file:
        try:
            panel_gene_dict[gene_type] = pd.read_csv(gene_file, sep="\t", header=None, comment="#").squeeze("columns")
        except pd.errors.EmptyDataError:
            panel_gene_dict[gene_type] = pd.Series(dtype='str')
    else:
        panel_gene_dict[gene_type] = None


fusion_table = pd.read_csv(args.input, sep="\t", header=0)
initial_column_list = fusion_table.columns

if panel_gene_dict["target"] is None:
    if panel_gene_dict["control"] is None:
        # No target or control genes were set
        fusion_table["is_target"] = True
        fusion_table["is_control"] = False
        fusion_table["is_offtarget"] = False
    else:
        fusion_table["is_control"] = fusion_table[args.gene1_column_name].isin(panel_gene_dict["control"]) | fusion_table[args.gene2_column_name].isin(panel_gene_dict["control"])
        fusion_table["is_target"] = ~fusion_table["is_control"]
        fusion_table["is_offtarget"] = False

else:
    if panel_gene_dict["control"] is None:
        fusion_table["is_control"] = False
        fusion_table["is_target"] = fusion_table[args.gene1_column_name].isin(panel_gene_dict["target"]) | fusion_table[args.gene2_column_name].isin(panel_gene_dict["target"])
        fusion_table["is_offtarget"] = ~fusion_table["is_target"]
    else:
        fusion_table["is_target"] = fusion_table[args.gene1_column_name].isin(panel_gene_dict["target"]) | fusion_table[args.gene2_column_name].isin(panel_gene_dict["target"])
        fusion_table["is_control"] = fusion_table[args.gene1_column_name].isin(panel_gene_dict["control"]) | fusion_table[args.gene2_column_name].isin(panel_gene_dict["control"]) & (~fusion_table["is_target"])
        fusion_table["is_offtarget"] = (~fusion_table["is_target"]) & (~fusion_table["is_control"])

fusion_table["gene_type"] = ""

for gene_type in "target", "control", "offtarget":
    fusion_table.loc[fusion_table["is_" + gene_type], "gene_type"] = gene_type

fusion_table[list(initial_column_list) + ["gene_type"]].to_csv(args.all_output, sep="\t", index=False, header=True)

fusion_table.loc[fusion_table["is_target"], initial_column_list].to_csv(args.target_output, sep="\t",
                                                                        index=False, header=True)
fusion_table.loc[fusion_table["is_control"], initial_column_list].to_csv(args.control_output, sep="\t",
                                                                         index=False, header=True)
fusion_table.loc[fusion_table["is_offtarget"], initial_column_list].to_csv(args.offtarget_output, sep="\t",
                                                                           index=False, header=True)






