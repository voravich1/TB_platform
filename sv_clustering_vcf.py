#!/usr/bin/env python

"""This program is use for cluster SV event form multiple vcf (like cluster program from genomadSV)
    Use exact mapping criteria
    May move to new project name SV_platform
"""

import sys
import json
import vcfpy
from collections import OrderedDict
import argparse
import os
import typing
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

mapping_mode = "exact"      # mapping criteria mode. Have 2 mode "exact" and "overlap"
sv_event_collate_dict = typing.Dict[str, list()]
sample_record_dict = dict()

for filename in os.listdir(input_directory):
    if filename.endswith(".vcf") or filename.endswith(".vcf.gz"):
        vcf_file = os.path.join(input_directory,filename)
        reader = vcfpy.Reader.from_path(vcf_file)

        samplename = vcf_file.split('_')[0].split('.')[0]   # two step split in case that sample name was separate by _ and .

        for record in reader:

            if(mapping_mode == "exact"):

                if(record.FILTER.contain("PASS") and "IMPRECISE" not in record.INFO):   # Filter only trust sv event
                    start = record.POS
                    info = record.INFO
                    end = info['END']
                    sv_type = info['SVTYPE']
                    sv_event_id = sv_type + "_" + str(start) + "_" + str(end)
                    sample_record_dict[samplename] = record

                    if(sv_event_id in sv_event_collate_dict):
                        sample_record_list = list()
                        sample_record_list = sv_event_collate_dict[sv_event_id]
                        sample_record_list
                    sv_event_collate_dict[sv_event_id] = sample_record_dict


