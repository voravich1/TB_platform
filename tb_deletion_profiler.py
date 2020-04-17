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
import ntpath
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

def binary_search_sv_mapping(in_vcf_record,in_del_db_dict_bisect, percent_overlap):
    # extract data from vcf record
    id = in_vcf_record.ID
    chr = "1" #in_vcf_record.CHROM
    info = in_vcf_record.INFO
    sv_type = info["SVTYPE"]
    result = dict()
    if "IMPRECISE" in info:  # cut event that has IMPRECISE flag
        return None
    elif sv_type != "DEL":
        return None

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
        ##index_end = bisect.bisect_left(list_sorted_adjust_overlap_end, end)
        ##dummy_end_info_pass_list = list_info_end[index_end:]

        # check start position. get all index that are lower or equal to query number
        # special treat on index on the left. we check equality of number index with query number
        # If it equal we will count it in. So, we plus 1 to index beacuase when we query from list with rage index it will no inclusive for back index. eg. num[:2] you will get value on 0 and 1 index not 0,1 and 2
        ##index_start = bisect.bisect_left(list_sorted_adjust_overlap_start, start)
        ##dummy_start_info_pass_list = list()
        ##if index_start == 12:
            ##print(index_start)

        #if index_start >= len(list_sorted_adjust_overlap_start): # in case that query start has value higher than last member bisect will return index that exceed the actual size of list which give error
                #dummy_start_info_pass_list = list_info_start
       #else:
            #if list_sorted_adjust_overlap_start[index_start] == start:
                #if len(list_info_start) != index_start:
                    #dummy_start_info_pass_list = list_info_start[:index_start+1]
                #else:
                    #dummy_start_info_pass_list = list_info_start
            #else:
                #dummy_start_info_pass_list = list_info_start[:index_start]

        # check end position. get all index that are higher than or equal to query number
        index_start = bisect.bisect_left(list_sorted_adjust_overlap_start, start)
        dummy_start_info_pass_list = list_info_start[index_start:]

        # check start position. get all index that are lower or equal to query number
        # special treat on index on the left. we check equality of number index with query number
        # If it equal we will count it in. So, we plus 1 to index beacuase when we query from list with rage index it will no inclusive for back index. eg. num[:2] you will get value on 0 and 1 index not 0,1 and 2
        index_end = bisect.bisect_left(list_sorted_adjust_overlap_end, end)
        dummy_end_info_pass_list = list()
        #if start == 4092079:
            #print(index_end)

        if index_end >= len(list_sorted_adjust_overlap_end):  # in case that query start has value higher than last member bisect will return index that exceed the actual size of list which give error
            dummy_end_info_pass_list = list_info_end
        else:
            if list_sorted_adjust_overlap_end[index_end] == end:
                if len(list_info_end) != index_end:
                    dummy_end_info_pass_list = list_info_end[:index_end + 1]
                else:
                    dummy_end_info_pass_list = list_info_end
            else:
                dummy_end_info_pass_list = list_info_end[:index_end]

        ###############################
        # loop check RD match pairing
        pairing_dict = dict()

        # loop start pass list
        for info in dummy_start_info_pass_list:
            rd_id = info[0]
            lineage = info[1]
            sv_len = info[4]

            pairing_dict[rd_id] = info

        # loop end pass list (judgement loop)
        result_count = 0
        for info in dummy_end_info_pass_list:
            rd_id = info[0]

            if rd_id in pairing_dict:
                start2 = info[2]
                end2 = info[3]
                ## Final check. Reciprocal overlap
                check_flag = check_reciprocal_overlap(start,end,start2,end2,percent_overlap)

                if check_flag == True:
                    #result_count+=1
                    result[info[0]] = info
        if not result:
            return None
        else:
            return result
    else:
        return None
    ####################################
########################################################################################################################

########################################################################################################################
## Check reciprocal overlap function
def check_reciprocal_overlap(startA, endA, startB, endB, percent_overlap):


    #event1 = list(range(startA, (endA + 1)))
    #event2 = list(range(startB, (endB + 1)))

    #region_intersect = intersection(event1,event2)
    #num_intersect = len(region_intersect)

    #event1_percent_overlap = (num_intersect * 100) / len(event1)
    #event2_percent_overlap = (num_intersect * 100) / len(event2)

    ##############################################################
    ## Numeric calculate overlap

    len_A = (endA - startA) + 1
    len_B = (endB - startB) + 1
    if startA >= startB and endA > endB:
        x = startA - startB
        num_overlap = len_B - x
    elif startA < startB and endA <= endB:
        x = startB - startA
        num_overlap = len_A - x
    elif startA >= startB and endA <= endB:
        x = startA - startB
        y = endB - endA
        u = len_B - x
        num_overlap = u - y
    elif startA < startB and endA > endB:
        x = startB - startA
        y = endA - endB
        u = len_A - x
        num_overlap = u - y
    elif startA == startB and endA == endB:
        num_overlap = len_A

    event1_percent_overlap_num = (num_overlap*100)/len_A
    event2_percent_overlap_num = (num_overlap*100)/len_B
    #########################################################

    #if num_overlap != num_intersect:
        #print("Noooooo!!!!\n")
        #print("starA:" + startA + "\n")
        #print("endA:" + endA + "\n")
        #print("starB:" + startB + "\n")
        #print("endB:" + endB + "\n")


    if event1_percent_overlap_num >= percent_overlap and event2_percent_overlap_num >= percent_overlap:
        return True
    else:
        return False
########################################################################################################################

########################################################################################################################
## Intersect two list
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3
########################################################################################################################

########################################################################################################################
## Phase one lineage classify decission tree (Hard code)
def phase_one_lineage_judgement(result_dict):
    lineage_result = dict()
    if "RD149" in result_dict:
        if "RD750" in result_dict:
            info1 = result_dict["RC149"]
            info2 = result_dict["RD750"]
            info = [info1,info2]
            lineage_result["lineage3"] = info
        elif "RD105_EX" in result_dict:
            info1 = result_dict["RC149"]
            info2 = result_dict["RD105_EX"]
            info = [info1, info2]
            lineage_result["lineage2-proto"] = info
        elif "RD105" in result_dict:
            info1 = result_dict["RC149"]
            info2 = result_dict["RD105"]
            info = [info1, info2]
            lineage_result["lineage2"] = info
        else:
            info = result_dict["RD149"]
            lineage_result["lineage4"] = info
    else:
        if "RD239" in result_dict:
            info = result_dict["RD239"]
            lineage_result["lineage1"] = info
            if "RV0209" in result_dict and "RV1004" in result_dict and "RV2531" in result_dict:
                info1 = result_dict["RV0209"]
                info2 = result_dict["RV01004"]
                info3 = result_dict["RV2531"]
                info = [info1,info2,info3]
                lineage_result["lineage 1.2.1"] = info
                if "RD121" in result_dict:
                    info = result_dict["RD121"]
                    lineage_result["lineage1.2.1.1"] = info
                elif "RV968" in result_dict:
                    info = result_dict["RV968"]
                    lineage_result["lineage 1.2.1.2"] = info
        elif "RD711" in result_dict:
            info = result_dict["RD711"]
            lineage_result["lineage5"] = info
        elif "RD702" in result_dict:
            info = result_dict["RD702"]
            lineage_result["lineage6"] = info
        else:
            lineage_result["unclassify_lineage"] = False

    return lineage_result
########################################################################################################################

vcfFile_name = args.input
delDBFile_name = args.lineage_del_db

reader_sv_vcf = vcfpy.Reader.from_path(vcfFile_name)

vcfFile_basename = ntpath.basename(vcfFile_name)
samplename = vcfFile_basename.split("_")[0]
#samplename = reader_sv_vcf.header.samples.names[0] ## this will return list of sample name in a header row of vcf

output_json_file_name = samplename + "_result_test.json"
json_result_file = os.path.join(args.output, output_json_file_name)

header_vcf = reader_sv_vcf.header
percent_overlap = 70

#del_db_dict = import_del_db(delDBFile_name)
del_db_bisect_dict = import_del_db_bisect(delDBFile_name,percent_overlap)

result_dict = dict()
result_dict_json = dict()
#result_count = 0
for record in reader_sv_vcf :
    #dummy_res_dict = sv_database_mapping(record,del_db_dict)
    dummy_res_dict = binary_search_sv_mapping(record, del_db_bisect_dict, percent_overlap)

    if dummy_res_dict != None:
        #result_dict +=1
        #result_dict[result_dict] = dummy_res_dict
        result_dict.update(dummy_res_dict)

lineage_final_result_dict = phase_one_lineage_judgement(result_dict)
lineage_final_result_sorted_dict = OrderedDict(sorted(lineage_final_result_dict.items()))

result_dict_json["sample_name"] = samplename
result_dict_json["lineage"] = lineage_final_result_sorted_dict

js = json.dumps(result_dict_json, sort_keys=True, indent=4)


json_write = open(json_result_file,'w')

json_write.write(js)

json_write.close()
print("Done read vcf")

    ## Test with single sample for RD239 pass but still has problem with false call lineage 4 that cause false call lineage1 insample ERR718222








