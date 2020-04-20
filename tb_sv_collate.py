#!/usr/bin/env python

"""This program is use for collate all json file result from deletion_profiler
    Compatible with json file from deletion_profiler only
"""

import sys
import os
import sys
import argparse
import json
import vcfpy
from collections import OrderedDict
from pathlib import Path

__author__ = "Worawich Phornsiricharoenphant"
__copyright__ = ""
__credits__ = ["Worawich Phornsiricharonphant"]
__license__ = "GPL-3.0"
__version__ = "1.0"
__maintainer__ = "Worawich Phornsiricharoenphant"
__email__ = "worawich.ph@gmail.com"
__status__ = "Development"


parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i', '--input', help='Input folder name', required=True)
parser.add_argument('-o', '--output', help='Absolute path of Output folder', required=True)
parser.add_argument('-t', '--target', help='Extension of target input file Ex. ".json"', default=".json", required=False)

args = parser.parse_args()


def collate_result(inputs, output, file_target=args.target):
    list_file = get_in_file_list(inputs, file_target)
    for file in list_file:
        with open(file) as fp:

            data_dict = json.load(fp)

            sample_name = data_dict["sample_name"]
            lineage_final_result = data_dict["lineage"]
            candidate_result = data_dict["candidate_results"]

#### continue create collate result of both snp and del in tab format for excel open