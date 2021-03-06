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
#drugDB = os.path.join(this_file_dir,"db","tbprofiler_drugDB.json")
drugDB = os.path.join(this_file_dir,"db","thai_drug_db.json")
linDB = os.path.join(this_file_dir,"db","lin_db_3916.txt")
#linDB = os.path.join(this_file_dir,"db","diff_maf_db")
#linDB = os.path.join(this_file_dir,"db","lin_tb_profiler.csv")

#vcfFile_snpeff_name = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/test/98_typing_snp_snpeff.vcf"      # vcf must be snnotate with snpEff

#drugDBFile_name = "/Users/worawich/Download_dataset/TB_platform_test/drug_db/tbprofiler_drugDB.json"        # need file in json. format inside will be follow tbprfiler drug db format

#lineageDBFile_name = "/Volumes/10TBSeagateBackupPlus/NBT/TB_platform/database/lineage_db/lin_db_6915.txt"       # need file in tab dilimit format [pos|ref|alt|lineage] no header

#json_result_file = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/test/98_lineage_drug_resultV4.json"     # result will be save in json format

#lineage_decision_threshold = 0.9        # lineage will be ofiicialy call as hit when 90% of marker was hit on each lineage
lineage_decision_threshold_mode = True     # True mean use threshold to decide lineage
lineage_decision_majority_mode = False      # True mean use majority vote to decide lineage
force_delete_non_majority_lineage = True    # True mean. At the last step final lineage that not the same as majority lineage will be delete

parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i', '--input', help='Input VCF file name (can be vcf.gz)', required=True)
parser.add_argument('-o', '--output', help='Absolute path of Output folder', required=True)
parser.add_argument('-d', '--drdb', help='adsolute path to drug database file (in .json format like tbprofiler used)', nargs='?', const=1, default=drugDB, required=False)
parser.add_argument('-l', '--lindb', help='adsolute path to lineage database file (in tab dilimit format [pos\tref\talt\tlineage])', nargs='?', const=1, default=linDB, required=False)
parser.add_argument('-e', '--threshold', help='lineage decision threshold value 0 to 1 default 0.9 which refer to 90%', const=1, nargs='?', type=float, default=0.9, required=False)

args = parser.parse_args()

vcfFile_snpeff_name = args.input
drugDBFile_name = args.drdb
lineageDBFile_name = args.lindb
lineage_decision_threshold = args.threshold    # lineage will be ofiicialy call as hit when 90% of marker was hit on each lineage

reader_snpeff_vcf = vcfpy.Reader.from_path(vcfFile_snpeff_name)

samplename = reader_snpeff_vcf.header.samples.names[0] ## this will return list of sample name in a header row of vcf

output_json_file_name = samplename + "_result_test.json"
json_result_file = os.path.join(args.output, output_json_file_name)

header_snpEff_vcf = reader_snpeff_vcf.header
snpEff_header = header_snpEff_vcf.get_info_field_info('ANN')
snpEff_header_description = snpEff_header.description
snpEff_field_order = snpEff_header_description.split(': \'')[1].split(' | ')
snpEff_field_query_index = dict()
hgvsc_snpEff_idx = 0
hgvsp_snpEff_idx = 0
hgvsc_snpEff_dict = dict()
hgvsp_snpEff_dict = dict()

annotation_snpEff_idx = 0

drugDB_gene_dict = dict()
drugDB_hgvs_dict = dict()
drugDB_drug_dict = dict()

special_drugDB_gene_dict = dict()
special_drugDB_hgvs_dict = dict()
special_drugDB_drug_dict = dict()

drug_result_dict = dict()

lineage_result_dict = dict()
lineage_final_result_dict = dict()
lineage_final_result_sorted_dict = dict()

result_dict = dict()

###########
## Loop check hgvs field index in snpEff vcf
###########
count = 0
for snpEff_field in snpEff_field_order:
    snpEff_field_query_index[snpEff_field] = count
    if(snpEff_field == 'HGVS.c'):
        hgvsc_snpEff_idx = count
    elif(snpEff_field == 'HGVS.p'):
        hgvsp_snpEff_idx = count
    elif(snpEff_field == 'Annotation'):
        annotation_snpEff_idx = count
    count+=1

#################################################

##########
## Loop file drug database
## drug database of TBprofiler is json format
## transform DB structure to NBT_AON structure [Gene => HGVS => Drug => Confident]
##########
aon_drug_list = list()
with open(drugDBFile_name) as json_file:
    data = json.load(json_file)

    for gene in data:

        gene_dict = data[gene]

        drugDB_hgvs_dict = dict()
        special_drugDB_hgvs_dict = dict()
        specialcase_flag = False
        for variant in gene_dict:

            #if variant == '-15C>T':
                #print("")

            variant_dict = gene_dict[variant]

            drug_dict = variant_dict['drugs']
            hgvs = variant_dict['hgvs_mutation']

            #if hgvs == "p.Ser450Leu": #"p.Lys43Arg" p.Leu452Pro
                #print()
            ## cheek special case that did not use hgvs format
            hgvs_split = hgvs.split(".")
            if(len(hgvs_split) == 1):
                # hgvs variable is not hgvs format. It's a special case
                special_drugDB_drug_dict = dict()
                for drug_name in drug_dict:
                    drug_name_dict = drug_dict[drug_name]
                    confidence = drug_name_dict['confidence']
                    if "resistance level" in drug_name_dict:
                        resistant_level = drug_name_dict['resistance level']
                    else:
                        resistant_level = "unknown"
                    dummy_dict = {"resistance_level":resistant_level,"confidence":confidence}
                    special_drugDB_drug_dict[drug_name] = dummy_dict

                    if drug_name not in aon_drug_list:
                        aon_drug_list.append(drug_name)

                special_drugDB_hgvs_dict[hgvs] = special_drugDB_drug_dict
                specialcase_flag = True
            #########################################################
            ## Normal case
            else:
                drugDB_drug_dict = dict()
                for drug_name in drug_dict:

                    drug_name_dict = drug_dict[drug_name]
                    confidence = drug_name_dict['confidence']
                    if "resistance level" in drug_name_dict:
                        resistant_level = drug_name_dict['resistance level']
                    else:
                        resistant_level = "unknown"
                    dummy_dict = {"resistance_level": resistant_level, "confidence": confidence}
                    drugDB_drug_dict[drug_name] = dummy_dict

                    if drug_name not in aon_drug_list:
                        aon_drug_list.append(drug_name)

                drugDB_hgvs_dict[hgvs] = drugDB_drug_dict
            ##############################################################
        drugDB_gene_dict[gene] = drugDB_hgvs_dict

        if(specialcase_flag == True): ## put special case info to special case dict
            special_drugDB_gene_dict[gene] = special_drugDB_hgvs_dict

json_file.close()
###########################################################

##########
## Loop file lineage database
## lineage data file should be in tap delimit format [POS\tREF\tALT\tLIN]
## No header in file
## Transform info to dict for later checking
##########
lin_db_dict = dict()
lin_db_total_dict = dict()
with open(lineageDBFile_name,'r') as file:
    for line in file:
        data = line.splitlines()[0].split('\t')
        pos = data[0]
        ref = data[1]
        alt = data[2]
        lin = data[3]
        key = pos + ref + alt
        lin_db_dict[key] = lin

        if lin in lin_db_total_dict:  ## add hit count to specific lineage
            count = lin_db_total_dict[lin]
            count += 1
            lin_db_total_dict[lin] = count
        else:
            lin_db_total_dict[lin] = 1
file.close()
###########################################################

###########
## For snpEff_VCF
## Loop Filter and Collect only record that give HGVS field
## Goal for loop on snpEff VCF is to get drug information
## (My understanding right now is "record that has HGVSc and HGVSp is variant record that occur on gene")
###########

record_hit_ID = 0
lin_hit_dict = dict()

for record in reader_snpeff_vcf:

    #########################################
    #### Lineage classify part
    #########################################

    pos = str(record.POS)
    ref_allel = record.REF
    alt_list = record.ALT

    if(len(alt_list) == 1): ## skip variant that have multiallele
        alt_allel = alt_list[0].value
        check_key = pos + ref_allel + alt_allel

        if check_key in lin_db_dict:    ## check variant hit on lineage DB
            lin = lin_db_dict[check_key]

            if lin in lin_hit_dict:     ## add hit count to specific lineage
                count = lin_hit_dict[lin]
                count+=1
                lin_hit_dict[lin] = count
            else:
                lin_hit_dict[lin] = 1

    #################################################################

    #########################################
    #### Drug resistant classify Part
    #########################################

    info = record.INFO
    record_ann = info.get('ANN')


    for snpEff_field in record_ann:
        snpEff_field_list = snpEff_field.split("|")

        if(snpEff_field_list[hgvsc_snpEff_idx] == "" and snpEff_field_list[hgvsp_snpEff_idx] == ""):
            continue

        snpEff_geneID = snpEff_field_list[snpEff_field_query_index['Gene_ID']]
        snpEff_geneName = snpEff_field_list[snpEff_field_query_index['Gene_Name']]
        annotation_snpEff = snpEff_field_list[annotation_snpEff_idx]

        ## Because we know how special case look like. So,reformat infomation to map with special case (hard code)
        ## there is 3 special case right now I just hard code on 2 case. the The third case wich is large_deletion sitll not be handle yet.
        special_case_check = "no_special_case"
        if(annotation_snpEff == "missense_variant"):
            dummy_hgvsp = snpEff_field_list[hgvsp_snpEff_idx]
            codon_number = snpEff_field_list[13].split("/")[0]
            if dummy_hgvsp != "":
             #   dummy_ref = record.REF
              #  split_word = dummy_ref + ">"
               # dummy_hgvsp_split = dummy_hgvsp.split(split_word)[0]
                #if len(dummy_hgvsp.split(split_word)) == 1: ## in case that hgvs is use complement base of reference
                 #   if dummy_ref == "A":
                  #      dummy_ref = "T"
                   # elif dummy_ref == "T":
                    #    dummy_ref = "A"
                    #elif dummy_ref == "C":
                    #    dummy_ref = "G"
                    #elif dummy_ref == "G":
                    #    dummy_ref = "C"
                    #split_word = dummy_ref + ">"
                    #dummy_hgvsp_split = dummy_hgvsp.split(split_word)[0]
                    #if dummy_hgvsp.split(split_word) == 1:
                        #print("")

                #dummy_nucleotide_number = dummy_hgvsp_split.split(".")[1]
                special_case_check = "any_missense_codon_" + codon_number
            else:
                special_case_check = "any_missense"
        elif(annotation_snpEff == "frameshift_variant"):
            special_case_check = "frameshift"
        ##

        hgvsc_info = snpEff_field_list[hgvsc_snpEff_idx]
        hgvsp_info = snpEff_field_list[hgvsp_snpEff_idx]
        ### Reformat hgvsc into the same format used in Database
        if hgvsp_info == "":    ## case check for rrs and rrl gene
            hgvsc_ready = hgvsc_info
            if snpEff_geneName == "rrs" or snpEff_geneName == "rrl": ## force change hgvs format of n. to r. To make it mapab;e to db of this 2 gene
                hgvsc_ready = hgvsc_info.replace("n","r",1).lower() ## Force transform lower case because base in data base is lower case
        elif len(hgvsc_info.split("del")) > 1: ## case check for small deletion events
            dummy_hgvsc_split=hgvsc_info.split("del")
            if dummy_hgvsc_split[1] != "":
                dummy_unuse_part=dummy_hgvsc_split[1]
                hgvsc_ready=hgvsc_info.split(dummy_unuse_part)[0]
            else:
                hgvsc_ready=hgvsc_info

        else:
            hgvsc_ready = hgvsc_info

        #########################################################
        if hgvsp_info != "":
            hgvsp_ready = hgvsp_info
        else:
            hgvsp_ready = hgvsp_info

        if hgvsp_ready == "p.Ser450Leu": #"p.Lys43Arg" "p.Leu452Pro"
            print("")

        if hgvsp_ready == "p.Asp94Ala": #"p.Lys43Arg" "p.Leu452Pro"
            print("")

        if hgvsc_ready == 'c.1290delC':
            print(snpEff_geneID)

        if hgvsc_ready == 'c.-15C>T':
            print(snpEff_geneID)

        if hgvsp_ready == 'c.-15C>T':
            print(snpEff_geneID)

        #########################################
        #### Check with database
        #########################################
        ## Special case check not hgvs as key (need special hard code)
        ## Check only for "any_missense_codon_" special case
        if snpEff_geneID in special_drugDB_gene_dict:
            gene_hit_dict = special_drugDB_gene_dict[snpEff_geneID]

            record_hit_dict = dict()
            if special_case_check in gene_hit_dict and special_case_check != "frameshift": ## "frameshift" will not check here but will be check after normal case not hit
                record_hit_dict['Drug'] = gene_hit_dict[special_case_check]
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(snpEff_field_order)):
                    record_hit_dict[snpEff_field_order[i]] = snpEff_field_list[i]

                drug_result_dict[record_hit_ID] = record_hit_dict
                record_hit_ID += 1

        elif snpEff_geneName in special_drugDB_gene_dict:
            gene_hit_dict = special_drugDB_gene_dict[snpEff_geneName]

            record_hit_dict = dict()
            if special_case_check in gene_hit_dict and special_case_check != "frameshift": ## "frameshift" will not check here but will be check after normal case not hit
                record_hit_dict['Drug'] = gene_hit_dict[special_case_check]
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(snpEff_field_order)):
                    record_hit_dict[snpEff_field_order[i]] = snpEff_field_list[i]

                drug_result_dict[record_hit_ID] = record_hit_dict
                record_hit_ID += 1
        ################################################################################
        ## Normal case check hgvs as key
        if snpEff_geneID in drugDB_gene_dict:
            gene_hit_dict = drugDB_gene_dict[snpEff_geneID]

            record_hit_dict = dict()
            if hgvsc_ready in gene_hit_dict:

                record_hit_dict['Drug'] = gene_hit_dict[hgvsc_ready]
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(snpEff_field_order)):
                    record_hit_dict[snpEff_field_order[i]] = snpEff_field_list[i]

                drug_result_dict[record_hit_ID] = record_hit_dict
                record_hit_ID += 1

            elif hgvsp_ready in gene_hit_dict:

                record_hit_dict['Drug'] = gene_hit_dict[hgvsp_ready]
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(snpEff_field_order)):
                    record_hit_dict[snpEff_field_order[i]] = snpEff_field_list[i]

                drug_result_dict[record_hit_ID] = record_hit_dict
                record_hit_ID += 1
            else:
                ## Normal case not hit. So, check for "frameshift" special case
                if snpEff_geneID in special_drugDB_gene_dict:
                    gene_hit_dict = special_drugDB_gene_dict[snpEff_geneID]

                    record_hit_dict = dict()
                    if special_case_check in gene_hit_dict and special_case_check == "frameshift":  ## "frameshift" will not check here but will be check after normal case not hit
                        record_hit_dict['Drug'] = gene_hit_dict[special_case_check]
                        record_hit_dict['Ref'] = record.REF
                        record_hit_dict['Pos'] = record.POS

                        for i in range(len(snpEff_field_order)):
                            record_hit_dict[snpEff_field_order[i]] = snpEff_field_list[i]

                        drug_result_dict[record_hit_ID] = record_hit_dict
                        record_hit_ID += 1
                #######################################################
        elif snpEff_geneName in drugDB_gene_dict:
            gene_hit_dict = drugDB_gene_dict[snpEff_geneName]

            record_hit_dict = dict()
            if hgvsc_ready in gene_hit_dict:

                record_hit_dict['Drug'] = gene_hit_dict[hgvsc_ready]
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(snpEff_field_order)):
                    record_hit_dict[snpEff_field_order[i]] = snpEff_field_list[i]

                drug_result_dict[record_hit_ID] = record_hit_dict
                record_hit_ID += 1

            elif hgvsp_ready in gene_hit_dict:

                record_hit_dict['Drug'] = gene_hit_dict[hgvsp_ready]
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(snpEff_field_order)):
                    record_hit_dict[snpEff_field_order[i]] = snpEff_field_list[i]

                drug_result_dict[record_hit_ID] = record_hit_dict
                record_hit_ID += 1
            else:
                ## Normal case not hit. So, check for "frameshift" special case
                if snpEff_geneName in special_drugDB_gene_dict:
                    gene_hit_dict = special_drugDB_gene_dict[snpEff_geneName]

                    record_hit_dict = dict()
                    if special_case_check in gene_hit_dict and special_case_check == "frameshift":  ## "frameshift" will not check here but will be check after normal case not hit
                        record_hit_dict['Drug'] = gene_hit_dict[special_case_check]
                        record_hit_dict['Ref'] = record.REF
                        record_hit_dict['Pos'] = record.POS

                        for i in range(len(snpEff_field_order)):
                            record_hit_dict[snpEff_field_order[i]] = snpEff_field_list[i]

                        drug_result_dict[record_hit_ID] = record_hit_dict
                        record_hit_ID += 1
                #########################################################
        #################################################################

        if hgvsc_ready in hgvsc_snpEff_dict:
            #print("repeat!!" + hgvsc_ready)
            record_list = hgvsc_snpEff_dict[hgvsc_ready]
            record_list.append(record)
            hgvsc_snpEff_dict[hgvsc_ready] = record_list
        else:
            record_list = [record]
            hgvsc_snpEff_dict[hgvsc_ready] = record_list

        if hgvsp_ready in hgvsp_snpEff_dict:
            #print("repeat!!" + hgvsp_ready)
            record_list = hgvsp_snpEff_dict[hgvsp_ready]
            record_list.append(record)
            hgvsp_snpEff_dict[hgvsp_ready] = record_list
        else:
            record_list = [record]
            hgvsp_snpEff_dict[hgvsp_ready] = record_list
    #########################################################
#########################################################

print("Done read vcf")

#########################################################

######################################
## Decide which lineage belong to this sample
######################################

if lineage_decision_threshold_mode == True:
    # Iterate over all the items in dictionary to find keys with max value
    # key is lineage, value is hit count

    # this flag use for check weather lineage4 is in candidate or not.
    # True mean => force assign lineage4
    # False mean => no need to force assign lineage4 because it already in candidate. Wether it pass or fail the judgment we are not right to force asiign it after that.
    force_assign_lineage4 = True
    force_assign_lineage4_9 = True
    ################################################################

    for key, value in lin_hit_dict.items():
        total_marker = lin_db_total_dict[key]

        ## check hit count with threshold to decide final lineage result
        if key == "lineage4":   ## some special lineage need special treat to re calculate hit count [lineage4]
            force_assign_lineage4 = False # Turn flag false => no need to assign lineage4 after this
            special_count = total_marker - value
            if special_count >= (total_marker * lineage_decision_threshold):
                lineage_final_result_dict[key] = special_count
        elif key == "lineage4.9":     ## some special lineage need special treat to re calculate hit count [lineage4.9] not compatible with tbprofiler DB
            force_assign_lineage4_9 = False  # Turn flag false => no need to assign lineage4.9 after this
            special_count = total_marker - value
            if special_count >= (total_marker * lineage_decision_threshold):
                lineage_final_result_dict[key] = special_count
        else: ## normal lineage
            if value >= (total_marker * lineage_decision_threshold):
                lineage_final_result_dict[key] = value

    ## Force assign lineage4 process
    if "lineage4" not in lin_db_total_dict:
        # last case check if no lineage4 in DB we will skip force assign process
        force_assign_lineage4 = False

    if force_assign_lineage4 == True:
        # Check for other main lineage and sub lineage hit except lineage4 family. If there is any main lineage or sub lineage already we will not force assign lineage 4 to result
        # the status is indicate by non_lineage4_hit_flag and non_lineag4_candidate_flag. [Hard code]
        non_lineage4_hit_flag = False
        lineage4_hit_flag = False
        non_main_lineage4_check_list = ["lineage1","lineage2","lineage3","lineage5","lineage6","lineage7"]
        for key in lineage_final_result_dict:
            dummy_lineage = key.split(".")[0]
            if dummy_lineage in non_main_lineage4_check_list:
                non_lineage4_hit_flag = True

        ## We focus on the case that has no main lineage or sub lineage assigned rather yet. Idicate by non_lineage_hit_flag.
        ## Force assign for lineage4 and lineage 4.9 (which is reference that we use now day, So in theory it should not match any snp marker that refer to lineage 4 in DB)
        ## quite a hard code, We give a full score which is a total number of marker have in DB for it. Becuase, normally we reverse the map score for this two lineage.
        ## As you can see it the case check above. [Noted: not sure lineage 4.9 should not be treat liek this or not]
        ## The another weak point of this hard fix is we can not detect mix sample with lineage4
        if non_lineage4_hit_flag == False:
            total_marker = lin_db_total_dict["lineage4"]
            lineage_final_result_dict["lineage4"] = total_marker
            #total_marker = lin_db_total_dict["lineage4.9"]
            #lineage_final_result_dict["lineage4.9"] = total_marker
    ####################################################################################################################

    ## Force assign lineage4.9 process
    if "lineage4.9" not in lin_db_total_dict:
        # last case check if no lineage4.9 in DB we will skip force assign process
        force_assign_lineage4_9 = False

    if force_assign_lineage4_9 == True:
        # Check for other main lineage and sub lineage hit except lineage4.9 family. If there is any main lineage or sub lineage already we will not force assign lineage 4.9 to result
        # the status is indicate by non_lineage4_9_hit_flag and non_lineag4_9_candidate_flag. [Hard code]
        non_lineage4_9_hit_flag = False
        lineage4_9_hit_flag = False
        non_lineage4_9_check_list = ["lineage1", "lineage2", "lineage3", "lineage5", "lineage6", "lineage7"]
        for key in lineage_final_result_dict:
            dummy_lineage = key.split(".")[0]
            if dummy_lineage in non_lineage4_9_check_list:
                non_lineage4_9_hit_flag = True

        ## We focus on the case that has no main lineage or sub lineage assigned rather yet. Idicate by non_lineage_hit_flag.
        ## Force assign for lineage4 and lineage 4.9 (which is reference that we use now day, So in theory it should not match any snp marker that refer to 4.9 in DB)
        ## quite a hard code, We give a full score which is a total number of marker have in DB for it. Becuase, normally we reverse the map score for this two lineage.
        ## As you can see it the case check above.
        ## The another weak point of this hard fix is we can not detect mix sample with lineage4.9
        if non_lineage4_9_hit_flag == False:
            total_marker = lin_db_total_dict["lineage4.9"]
            lineage_final_result_dict["lineage4.9"] = total_marker

    ####################################################################################################################


    lineage_final_result_sorted_dict = OrderedDict(sorted(lineage_final_result_dict.items())) ## sort dict by key
#########################################################

######################################
## Force delete non majority lineage result
## By scoring mainlineage hit give 10 score sub lineage hit in any level give 1 score
######################################
if len(lineage_final_result_sorted_dict) == 0:  ## In case that No lineage result assign we skip Force delete process
    force_delete_non_majority_lineage = False

if force_delete_non_majority_lineage == True:
    lineage_score_dict = dict()
    for key, value in lineage_final_result_sorted_dict.items():
        dummy_lineage_split = key.split(".")
        main_lineage = dummy_lineage_split[0]

        ## assign score
        score = 0
        if len(dummy_lineage_split) == 1:
            score = 10
        elif len(dummy_lineage_split) > 1:
            score = 1
        else:
            raise Exception('Has problem with Lineage String. It may empty?. Lineage string: {}'.format(key))
        ################

        ## put score to dict
        if main_lineage in lineage_score_dict:
            dummy_score = lineage_score_dict[main_lineage]
            update_score = dummy_score + score
            lineage_score_dict[main_lineage] = update_score
        else:
            lineage_score_dict[main_lineage] = score
        #####################

    ## Get highest score lineage list
    highest_item = max(lineage_score_dict.items(), key=lambda x: x[1])
    highest_score_lineage = highest_item[0]
    highest_score = highest_item[1]

    highest_score_lineage_list = list()
    for key, value in lineage_score_dict.items():
        if value == highest_score:
            highest_score_lineage_list.append(key)
    ##################################

    ## Remove non highest score from final result
    del_key_list = list()
    for key, value in lineage_final_result_sorted_dict.items():
        dummy_lineage_split = key.split(".")
        main_lineage = dummy_lineage_split[0]

        if main_lineage not in highest_score_lineage_list:
            del_key_list.append(key)
            #del lineage_final_result_sorted_dict[key]
    for key in del_key_list:
        del lineage_final_result_sorted_dict[key]
    ###############################################

######################################


######################################
## Select majority vote on lineage result
######################################

if lineage_decision_majority_mode == True:

    itemMaxValue = max(lin_hit_dict.items(), key=lambda x: x[1])

    listOfMaxKeys = list()
    # Iterate over all the items in dictionary to find keys with max value
    for key, value in lin_hit_dict.items():
        if value == itemMaxValue[1]:
            listOfMaxKeys.append(key)

    #print('Keys with maximum Value in Dictionary : ', listOfKeys)

    for lineage in lin_hit_dict:
        main_lineage = lineage.split(".")[0]
        if main_lineage in listOfMaxKeys:
           lineage_final_result_dict[lineage] = lin_hit_dict[lineage]

    lineage_final_result_sorted_dict = OrderedDict(sorted(lineage_final_result_dict.items()))
#########################################################

###########################################
## Classification of MDR and XDR (follow WHO)
###########################################

drug_all_hit = list()
for key in drug_result_dict:
    content = drug_result_dict[key]

    drug_hit = content['Drug']

    for drug in drug_hit:

        if drug not in drug_all_hit:
            drug_all_hit.append(drug)

if len(drug_all_hit) == 0:
    drug_resist_type = "sensitive"
else:
    if ('isoniazid' in drug_all_hit) and ("rifampicin" in drug_all_hit):

        # check fluoroquinolone derivative
        if ('ciprofloxacin' in drug_all_hit) or ('garenoxacin'  in drug_all_hit) or ('gatifloxacin' in drug_all_hit) or ('gemifloxacin' in drug_all_hit) or ('levofloxacin' in drug_all_hit) or ('moxifloxacin' in drug_all_hit):

            # check with secondline injectable drug
            if ('capreomycin' in drug_all_hit) or ('kanamycin' in drug_all_hit) or ('amikacin' in drug_all_hit):

                drug_resist_type = "XDR"
            else:
                drug_resist_type = "MDR"
        else:
            drug_resist_type = "MDR"
    else:
        drug_resist_type = "resistant"
###########################################

###########################################
## Create json dict contain lineage, drug result and drug resistant type
## Then save to json file
###########################################

result_dict["sample_name"] = samplename
result_dict["lineage"] = lineage_final_result_sorted_dict
result_dict["small_variant_dr"] = drug_result_dict
result_dict["drug_resist_type"] = drug_resist_type
js = json.dumps(result_dict, sort_keys=True, indent=4)


json_write = open(json_result_file,'w')

json_write.write(js)

json_write.close()
print("Done read vcf")


