#!/usr/bin/env python

"""This program is use for classify lineage fron SV event
    Have solid SV marker to classify main lineage but sub lineage still need more research on it.
    Compatible for vcf from manta or delly
"""

import sys
import json
import vcfpy
from collections import OrderedDict
from pathlib import Path

__author__ = "Worawich Phornsiricharoenphant,Wasna Viratyosin"
__copyright__ = ""
__credits__ = ["Worawich Phornsiricharonphant"]
__license__ = "GPL-3.0"
__version__ = "1.0"
__maintainer__ = "Worawich Phornsiricharoenphant"
__email__ = "worawich.ph@gmail.com"
__status__ = "Development"

sv_db_file = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/lin_db/sv_lin_db.txt"

with open(sv_db_file,'r') as sv_db_reader:
    for line in sv_db_reader:
