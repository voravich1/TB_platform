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
# Inner function for read deletion database file
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
# Inner method for mapping sv event to database
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

########################################################################################################################
# Inner function for read deletion database file convert to compatible for bisect mapping algorithm
def import_del_db_bisect(in_del_db_file, percent_overlap: int = 70):
    full_db_dict = dict()

    with open(in_del_db_file, 'r') as file:
        for line in file:
            data = line.splitlines()[0].split('\t')
            chr = data[0]
            start = int(data[1])
            end = int(data[2])
            sv_len = int(data[3])
            rd_id = data[4]
            lineage = data[5]
            predict_info = [rd_id, lineage, start, end, sv_len]


            adjust_pos = adjust_to_overlap_pos(start,end,sv_len,percent_overlap)
            adjust_start = adjust_pos[0]
            adjust_end = adjust_pos[1]

            if chr in full_db_dict:
                # get existing list
                dummy_main_db_dict = full_db_dict[chr]
                list_sorted_adjust_overlap_start = dummy_main_db_dict["START"]
                list_info_start = dummy_main_db_dict["START_INFO"]
                list_sorted_adjust_overlap_end = dummy_main_db_dict["END"]
                list_info_end = dummy_main_db_dict["END_INFO"]

                ## update list
                list_start_insert_index = bisect.bisect_left(list_sorted_adjust_overlap_start, adjust_start)
                list_sorted_adjust_overlap_start.insert(list_start_insert_index, adjust_start)
                list_info_start.insert(list_start_insert_index, predict_info)

                list_end_insert_index = bisect.bisect_left(list_sorted_adjust_overlap_end, adjust_end)
                list_sorted_adjust_overlap_end.insert(list_end_insert_index, adjust_end)
                list_info_end.insert(list_end_insert_index, predict_info)

                ## update dict
                dummy_main_db_dict["START"] = list_sorted_adjust_overlap_start
                dummy_main_db_dict["START_INFO"] = list_info_start
                dummy_main_db_dict["END"] =  list_sorted_adjust_overlap_end
                dummy_main_db_dict["END_INFO"] = list_info_end

                full_db_dict[chr] = dummy_main_db_dict
            else:
                # initial list dict and list
                dummy_main_db_dict = dict()
                list_sorted_adjust_overlap_start = [adjust_start]
                list_info_start = [predict_info]                    ## Has same order as list_sorted_adjust_overlap_start
                list_sorted_adjust_overlap_end = [adjust_end]
                list_info_end = [predict_info]                      ## Has same order as list_sorted_adjust_overlap_end

                dummy_main_db_dict["START"] = list_sorted_adjust_overlap_start
                dummy_main_db_dict["START_INFO"] = list_info_start
                dummy_main_db_dict["END"] = list_sorted_adjust_overlap_end
                dummy_main_db_dict["END_INFO"] = list_info_end

                full_db_dict[chr] = dummy_main_db_dict

    file.close()
    return full_db_dict
########################################################################################################################


########################################################################################################################
# Inner function for adjust start and end pos to accept percent overlap pos
def adjust_to_overlap_pos(in_start: int,in_end: int,event_size, percent_overlap):
    accept_ovelap_event_size = int((percent_overlap * event_size)/100)

    adjust_start = in_end - accept_ovelap_event_size
    adjust_end = in_start + accept_ovelap_event_size

    return [adjust_start,adjust_end]

########################################################################################################################


########################################################################################################################
# Inner method for mapping sv event to database (overlap logic)
import bisect

def binary_search_sv_mapping(in_vcf_record,in_del_db_dict_bisect):
    # extract data from vcf record
    id = in_vcf_record.ID
    chr = "1" #in_vcf_record.CHROM
    info = in_vcf_record.INFO
    result = dict()
    if "IMPRECISE" in info:  # cut event that has IMPRECISE flag
        return None

    sv_type = info["SVTYPE"]

    start = int(in_vcf_record.POS)  # 1 based coordinate the diff of 1 based and 0 based is for start 1 based is inclusive the number you saw but 0 based is not.
    end = int(info["END"])  # 1 based coordinate he diff of 1 based and 0 based is for For stop position both are inclusive
    if "SVLEN" in info:
        sv_len = int(info["SVLEN"][0])
    else:
        sv_len = "none"
    format = record.FORMAT
    sv_id = record.ID
    ####################################

    # extract data from db dict (must be bisect dict from import_del_db_bisect only)
    if chr in in_del_db_dict_bisect:
        main_db_dict = in_del_db_dict_bisect[chr]
        list_sorted_adjust_overlap_start = main_db_dict["START"]
        list_info_start = main_db_dict["START_INFO"]
        list_sorted_adjust_overlap_end = main_db_dict["END"]
        list_info_end = main_db_dict["END_INFO"]

        # bisect left will return insertion point or index on the left of match number. The number on that index are higher or equal than the query number
        # This mean, any index on the right (inclue the index it self) has value equal or higher than query number
        # and any index on the left has value lower than query number

        # check end position. get all index that are higher than or equal to query number
        index_end = bisect.bisect_left(list_sorted_adjust_overlap_end, end)
        dummy_end_info_pass_list = list_info_end[index_end:]

        # check start position. get all index that are lower or equal to query number
        # special treat on index on the left. we check equality of number index with query number
        # If it equal we will count it in. So, we plus 1 to index beacuase when we query from list with rage index it will no inclusive for back index. eg. num[:2] you will get value on 0 and 1 index not 0,1 and 2
        index_start = bisect.bisect_left(list_sorted_adjust_overlap_start, start)
        dummy_start_info_pass_list = list()
        if index_start == 12:
            print(index_start)

        if index_start >= len(list_sorted_adjust_overlap_start): # in case that query start has value higher than last member bisect will return index that exceed the actual size of list which give error
                dummy_start_info_pass_list = list_info_start
        else:
            if list_sorted_adjust_overlap_start[index_start] == start:
                if len(list_info_start) != index_start:
                    dummy_start_info_pass_list = list_info_start[:index_start+1]
                else:
                    dummy_start_info_pass_list = list_info_start
            else:
                dummy_start_info_pass_list = list_info_start[:index_start]
        ###############################
        # loop check RD match pairing
        pairing_dict = dict()

        # loop start pass list
        for info in dummy_start_info_pass_list:
            rd_id = info[0]
            lineage = info[1]
            start = info[2]
            end = info[3]
            sv_len = info[4]

            pairing_dict[rd_id] = info

        # loop end pass list (judgement loop)
        result_count = 0
        for info in dummy_end_info_pass_list:
            rd_id = info[0]

            if rd_id in pairing_dict:
                result_count+=1
                result[result_count] = info
        if not result:
            return None
        else:
            return result
    else:
        return None
    ####################################
########################################################################################################################



vcfFile_name = args.input
delDBFile_name = args.lineage_del_db

reader_sv_vcf = vcfpy.Reader.from_path(vcfFile_name)

samplename = reader_sv_vcf.header.samples.names[0] ## this will return list of sample name in a header row of vcf

output_json_file_name = samplename + "_result_test.json"
json_result_file = os.path.join(args.output, output_json_file_name)

header_vcf = reader_sv_vcf.header

#del_db_dict = import_del_db(delDBFile_name)
del_db_bisect_dict = import_del_db_bisect(delDBFile_name)

for record in reader_sv_vcf :
    #dummy_res_dict = sv_database_mapping(record,del_db_dict)
    dummy_res_dict = binary_search_sv_mapping(record, del_db_bisect_dict)
    print()
    ## continue code for finalize lineage need to wait for all record to be process







