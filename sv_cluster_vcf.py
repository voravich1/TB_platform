#!/usr/bin/env python

"""This program is use for cluster SV event form multiple vcf (like cluster program from genomadSV)
    Not finish
    move to new project name SV_platform
"""

import sys
import json
import vcfpy
from collections import OrderedDict
import argparse
import os
from pathlib import Path

__author__ = "Worawich Phornsiricharoenphant"
__copyright__ = ""
__credits__ = ["Worawich Phornsiricharonphant"]
__license__ = "GPL-3.0"
__version__ = "1.0"
__maintainer__ = "Worawich Phornsiricharoenphant"
__email__ = "worawich.ph@gmail.com"
__status__ = "Development"

input_directory = "/Volumes/10TBSeagateBackupPlus/NBT/sv_platform/data/data_test_sv_cluster"

for filename in os.listdir(input_directory):
    if filename.endswith(".vcf") or filename.endswith(".vcf.gz"):
        vcf_file = os.path.join(input_directory,filename)
        reader = vcfpy.Reader.from_path(vcf_file)

        for record in reader:
            if(record.FILTER.contain("PASS") and "IMPRECISE" not in record.INFO):
                start = record.POS
                info = record.INFO
                end = info['END']
