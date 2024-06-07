#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="File with fusions called by Arriba. Default: stdin")
parser.add_argument("-r", "--removed", action="store", dest="removed", default=None,
                    help="File with filtered out fusions. Default: not reported ")
parser.add_argument("-w", "--all_with_filters", action="store", dest="all_with_filters", default=None,
                    help="File with all fusions and filters shown. Default: not reported ")
parser.add_argument("-x", "--min_coverage", action="store", dest="min_coverage", default=1, type=int,
                    help="Minimal full coverage (coverage1 and coverage2) value to retain the fusion. Default: 1 ")
parser.add_argument("-n", "--min_non_zero_sup_cov", action="store", default=2, type=int,
                    dest="min_non_zero_sup_cov",
                    help="Minimal number of supporting read types (split_reads1,"
                         "split_reads2 and discordant_mates) above zero. "
                         "Allowed: 0, 1, 2(default), 3 ")
parser.add_argument("-a", "--coverage_ratio_threshold", action="store", dest="coverage_ratio_threshold",
                    default=0.001, type=float,
                    help="Threshold for coverage ratio filter. Default: 0.001 ")
parser.add_argument("-m", "--min_supporting_ratios_above_threshold", action="store",
                    dest="min_supporting_ratios_above_threshold",
                    default=2, type=int,
                    help="Minimal number of supporting ratios (split_reads1/coverage1,"
                         "split_reads2/coverage2 and max(discordant_mates/coverage)) above threshold. "
                         "Allowed: 0, 1, 2(default), 3 ")

parser.add_argument("-f", "--filtered", action="store", dest="filtered", default=sys.stdout,
                    help="File with filtered fusions. Default: stdout")
args = parser.parse_args()

if (args.min_supporting_ratios_above_threshold > 3) or (args.min_supporting_ratios_above_threshold < 0):
    raise ValueError(f"ERROR!!! Allowed values for --min_supporting_ratios_above_threshold are 0,1,2,3 only. "
                     " Value {args.min_supporting_ratios_above_threshold} is not allowed.")

fusion_table = pd.read_csv(args.input, sep="\t", header=0,)
initial_column_list = fusion_table.columns

# Checking if full coverage is zero in any of fused transcripts/genes
fusion_table["cov1_filter"] = fusion_table["coverage1"] >= args.min_coverage
fusion_table["cov2_filter"] = fusion_table["coverage2"] >= args.min_coverage

# Checking if at least two types
# of supporting reads ("split_reads1", "split_reads2" and "discordant_mates") are non zero
fusion_table["non_zero_sup_cov_filter"] = (fusion_table[["split_reads1", "split_reads2", "discordant_mates"]] > 0).sum(axis=1) >= args.min_non_zero_sup_cov

# Checking ration of supporting reads to the full coverages
fusion_table["(split_reads1/coverage1)_filter"] = (fusion_table["split_reads1"] / fusion_table["coverage1"]) >= args.coverage_ratio_threshold
fusion_table["(split_reads2/coverage2)_filter"] = (fusion_table["split_reads2"] / fusion_table["coverage2"]) >= args.coverage_ratio_threshold
fusion_table["max(discordant_mates/coverage)_filter"] = \
    (fusion_table["discordant_mates"]/ fusion_table["coverage1"]).compare(fusion_table["discordant_mates"] / fusion_table["coverage2"],
                                                                          keep_shape=True, keep_equal=True).max(axis=1) >= args.coverage_ratio_threshold

fusion_table["min_sup_cov_ratios_filter"] = fusion_table[["(split_reads1/coverage1)_filter",
                                                          "(split_reads2/coverage2)_filter",
                                                          "max(discordant_mates/coverage)_filter"]].sum(axis=1) >= args.min_supporting_ratios_above_threshold

filtered_fusion_series = fusion_table["cov1_filter"] & fusion_table["cov2_filter"] & \
                         fusion_table["non_zero_sup_cov_filter"] & fusion_table["min_sup_cov_ratios_filter"]
fusion_table["status"] = ""
fusion_table.loc[filtered_fusion_series, "status"] = "retained"
fusion_table.loc[~filtered_fusion_series, "status"] = "filtered_out"

filtered_fusion_table = fusion_table.loc[filtered_fusion_series, initial_column_list]
filtered_fusion_table.to_csv(args.filtered, sep="\t", header=True, index=False)

if args.all_with_filters is not None:
    all_with_filters_fusion_table = fusion_table
    all_with_filters_fusion_table.to_csv(args.all_with_filters, sep="\t", header=True, index=False)

if args.removed is not None:
    filtered_out_fusion_table = fusion_table.loc[~filtered_fusion_series, initial_column_list]
    filtered_out_fusion_table.to_csv(args.removed, sep="\t", header=True, index=False)