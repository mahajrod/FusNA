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


def find_files(directory, extension_list):
    dir_path = directory if isinstance(directory, PosixPath) else Path(directory)
    file_list = []
    for extension in [extension_list] if isinstance(extension_list, str) else extension_list:
        file_list += sorted(list(directory.glob("*{0}".format(extension))))
    return file_list

