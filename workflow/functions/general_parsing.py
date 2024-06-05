#!/usr/bin/env python
__author__ = "mahajrod"
"""
This file contains functions necessary for Snakemake file
"""
import os
import yaml
import shutil
from copy import deepcopy
from collections import OrderedDict
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path, PosixPath


def p_distance(seq_a, seq_b, seq_len):
    dist = 0
    for i in range(0, seq_len):
        if seq_a[i] != seq_b[i]:
            dist += 1
    return dist


def get_pair_prefix(forward_filename, reverse_filename):
    if len(forward_filename) != len(reverse_filename):
        raise ValueError(f"Paired files {forward_filename} and {reverse_filename} have different lengths of the filenames!")
    seq_len = len(forward_filename)
    end_prefix_index = 0
    for i in range(0, seq_len):
        if forward_filename[i] != reverse_filename[i]:
            end_prefix_index = i
            break
    if end_prefix_index == 0:
        raise ValueError(f"Paired files {forward_filename} and {reverse_filename} differ at first character!")
    return forward_filename[0:end_prefix_index]


def find_files(directory, extension_list):
    dir_path = directory if isinstance(directory, PosixPath) else Path(directory)
    file_list = []
    for extension in [extension_list] if isinstance(extension_list, str) else extension_list:
        file_list += sorted(list(dir_path.glob("*{0}".format(extension))))
    return file_list

"""
def find_files(directory, extension_list):
    dir_path = directory if isinstance(directory, PosixPath) else Path(directory)
    ext_list = [extension_list] if isinstance(extension_list, str) else extension_list
    return {extension: sorted(list(dir_path.glob("*{0}".format(extension)))) for extension in ext_list}
"""
