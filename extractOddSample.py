#!/usr/bin/env python

"""This program is use to extract odd sample from del extract file of other sv extract file
    loop to any SV extract file and pick sample that not have SV event which specific by target file

    input: 
    1. SV extract file format eg. del_extract.txt
    2. target SV evenet file in txt format 
        [single column SV event location the SV info must be the same ad appear in extract format eg. 1:1414787-1414948]
    
"""

import sys
import json
import argparse
import os
import re

__author__ = "Worawich Phornsiricharoenphant"
__copyright__ = ""
__credits__ = ["Worawich Phornsiricharonphant"]
__license__ = "GPL-3.0"
__version__ = "1.0"
__maintainer__ = "Worawich Phornsiricharoenphant"
__email__ = "worawich.ph@gmail.com"
__status__ = "Development"


input_extract_file="/Users/worawich/Downloads/TB_del_paper/del_analysis/sv_master_del_extract.txt"
input_sv_target_file=""
output_odd_file="/Users/worawich/Downloads/TB_del_paper/del_analysis/sv_master_del_extract_odd.txt"
header = True

header_list = []

odd_file = open(output_odd_file, "w")

## Read SV extract file
with open(input_extract_file,"r",encoding='utf8') as f:
    for line in f:
        if header == True:
            header_list = line.split("\t")
            header = False
            continue

        info_list = line.split("\t")
        count=0
        first_colume=True
        for info in info_list:
            if first_colume == True:
                odd_file.write(info)
                count+=1
                first_colume=False
                continue

            if info != "2":
                odd_sample = header_list[count]
                odd_file.write("\t")
                odd_file.write(odd_sample)
            count+=1    
        odd_file.write("\n")

odd_file.close()


