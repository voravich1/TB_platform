#!/usr/bin/env python

"""This program is use for classify lineage and drug resistant form vcf [single sample vcf] that already annotate by snpEff
    Suggest lineage DB is lin_db_6915 or lin_db_5654 (credit P'big mahidol, not yet ready to use)
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
__credits__ = ["Worawich Phornsiricharonphant","Wasna Viratyosin"]
__license__ = "GPL-3.0"
__version__ = "1.0"
__maintainer__ = "Worawich Phornsiricharoenphant"
__email__ = "worawich.ph@gmail.com"
__status__ = "Development"

this_file_dir = os.path.dirname(os.path.realpath(__file__))
del_db = os.path.join(this_file_dir,"db","tb_del_db_phase1.txt")

parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i', '--input', help='Input VCF file name (can be vcf.gz)', required=True)
parser.add_argument('-o', '--output', help='Absolute path of Output folder', required=True)
parser.add_argument('-d', '--lineage_del_db', help='adsolute path to lineage deletion marker database file (in tab delimit format)', nargs='?', const=1, default=del_db, required=False)

args = parser.parse_args()

########################################################################################################################
## Inner function for read deletion database file
def import_del_db(in_del_db_file):
    full_db_dict = dict()

    with open(in_del_db_file, 'r') as file:
        for line in file:
            data = line.splitlines()[0].split('\t')
            chr = data[0]
            start = data[1]
            end = data[2]
            sv_len = data[3]
            rd_id = data[4]
            lineage = data[5]

            predict_info = [rd_id, lineage, sv_len]

            if chr in full_db_dict:
                dummy_main_db_dict = full_db_dict[chr]
                if start in dummy_main_db_dict:
                    dummy_inner_db_dict = dummy_main_db_dict[start]
                    if end in dummy_inner_db_dict:
                        print("There is duplicate event in your Data base")
                    else:
                        predict_info = [rd_id,lineage,sv_len]
                        dummy_inner_db_dict[end] = predict_info
                    dummy_main_db_dict[start] = dummy_inner_db_dict
                else:
                    dummy_inner_db_dict = {end:predict_info}
                    dummy_main_db_dict[start] = dummy_inner_db_dict
                full_db_dict[chr] = dummy_main_db_dict
            else:
                dummy_inner_db_dict = {end:predict_info}
                dummy_main_db_dict = {start:dummy_inner_db_dict}
                full_db_dict = {chr:dummy_main_db_dict}
    file.close()
    return full_db_dict
########################################################################################################################


########################################################################################################################
## Inner method for mapping sv event to database
def sv_database_mapping(in_vcf_record,in_del_db_dict):        ## mapping sv event to DB. for first version we use exact match stratergy
    id = in_vcf_record.ID
    chr = in_vcf_record.CHROM
    info = in_vcf_record.INFO
    result = dict()
    if "IMPRECISE" in info:     # cut event that has IMPRECISE flag
        return None

    sv_type = info["SVTYPE"]

    start = in_vcf_record.POS  # 1 based coordinate the diff of 1 based and 0 based is for start 1 based is inclusive the number you saw but 0 based is not.
    end = info["END"]  # 1 based coordinate he diff of 1 based and 0 based is for For stop position both are inclusive
    if "SVLEN" in info:
        sv_len = info["SVLEN"]
    else:
        sv_len = "none"
    format = record.FORMAT
    sv_id = record.ID

    if chr in in_del_db_dict:
        main_del_db_dict = in_del_db_dict[chr]
        if start in main_del_db_dict:
            inner_del_db_dict = main_del_db_dict[start]
            if end in  inner_del_db_dict:
                rd_id = inner_del_db_dict[end]
                rd_id
                dummy_result = {"ID":id,"START":start,"END":end,"TYPE":sv_type, "LENGTH":sv_len}
                result[rd_id] = dummy_result
                return result
            else:
                return None
        else:
            return None
    else:
        return None

########################################################################################################################

vcfFile_name = args.input
delDBFile_name = args.lineage_del_db

reader_sv_vcf = vcfpy.Reader.from_path(vcfFile_name)

samplename = reader_sv_vcf.header.samples.names[0] ## this will return list of sample name in a header row of vcf

output_json_file_name = samplename + "_result_test.json"
json_result_file = os.path.join(args.output, output_json_file_name)

header_vcf = reader_sv_vcf.header

del_db_dict = import_del_db(delDBFile_name)

for record in reader_sv_vcf :
    dummy_res_dict = sv_database_mapping(record,del_db_dict)
    print()








