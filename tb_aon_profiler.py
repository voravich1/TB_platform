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
drugDB = os.path.join(this_file_dir,"db","tbprofiler_drugDB.json")
linDB = os.path.join(this_file_dir,"db","lin_db_6915.txt")

#vcfFile_snpeff_name = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/test/98_typing_snp_snpeff.vcf"      # vcf must be snnotate with snpEff

#drugDBFile_name = "/Users/worawich/Download_dataset/TB_platform_test/drug_db/tbprofiler_drugDB.json"        # need file in json. format inside will be follow tbprfiler drug db format

#lineageDBFile_name = "/Volumes/10TBSeagateBackupPlus/NBT/TB_platform/database/lineage_db/lin_db_6915.txt"       # need file in tab dilimit format [pos|ref|alt|lineage] no header

#json_result_file = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/test/98_lineage_drug_resultV4.json"     # result will be save in json format

#lineage_decision_threshold = 0.9        # lineage will be ofiicialy call as hit when 90% of marker was hit on each lineage
lineage_decision_threshold_mode = True     # True mean use threshold to decide lineage
lineage_decision_majority_mode = False      # True mea use majority vote to decide lineage

parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i', '--input', help='Input VCF file name (can be vcf.gz)', required=True)
parser.add_argument('-o', '--output', help='Absolute path of Output folder', required=True)
parser.add_argument('-d', '--drdb', help='adsolute path to drug database file (in .json format like tbprofiler used)', nargs='?', const=1, default=drugDB, required=False)
parser.add_argument('-l', '--lindb', help='adsolute path to lineage database file (in tab dilimit format [pos\tref\talt\tlineage])', nargs='?', const=1, default=linDB, required=False)
parser.add_argument('-e', '--threshold', help='lineage decision threshold value 0 to 1 default 0.9 which refer to 90%', const=1, nargs='?', type=int, default=0.9, required=False)

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

drugDB_gene_dict = dict()
drugDB_hgvs_dict = dict()
drugDB_drug_dict = dict()

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

    count+=1

#################################################

##########
## Loop file drug database
## drug database of TBprofiler is json format
## transform DB structure to NBT_AON structure [Gene => HGVS => Drug => Confident]
##########

with open(drugDBFile_name) as json_file:
    data = json.load(json_file)

    for gene in data:

        gene_dict = data[gene]

        drugDB_hgvs_dict = dict()
        for variant in gene_dict:

            if variant == '-15C>T':
                print("")

            variant_dict = gene_dict[variant]

            drug_dict = variant_dict['drugs']
            hgvs = variant_dict['hgvs_mutation']

            if hgvs == "p.Ser450Leu": #"p.Lys43Arg" p.Leu452Pro
                print()

            drugDB_drug_dict = dict()
            for drug_name in drug_dict:

                drug_name_dict = drug_dict[drug_name]
                confidence = drug_name_dict['confidence']

                drugDB_drug_dict[drug_name] = confidence

            drugDB_hgvs_dict[hgvs] = drugDB_drug_dict

        drugDB_gene_dict[gene] = drugDB_hgvs_dict
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

        hgvsc_info = snpEff_field_list[hgvsc_snpEff_idx]
        hgvsp_info = snpEff_field_list[hgvsp_snpEff_idx]

        if hgvsp_info != "":
            hgvsc_ready = hgvsc_info
        else:
            hgvsc_ready = hgvsc_info

        if hgvsp_info != "":
            hgvsp_ready = hgvsp_info
        else:
            hgvsp_ready = hgvsp_info

        if hgvsp_ready == "p.Ser450Leu": #"p.Lys43Arg" "p.Leu452Pro"
            print("")

        if hgvsc_ready == 'c.-15C>T':
            print(snpEff_geneID)

        if hgvsp_ready == 'c.-15C>T':
            print(snpEff_geneID)

        #########################################
        #### Check with database
        #########################################
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
    for key, value in lin_hit_dict.items():
        total_marker = lin_db_total_dict[key]

        ## check hit count with threshold to decide final lineage result
        if key == "lineage4" or key == "lineage4.9":   ## some special lineage need special treat to re calculate hit count
            special_count = total_marker - value
            if special_count >= (total_marker * lineage_decision_threshold):
                lineage_final_result_dict[key] = special_count
        else: ## normal lineage
            if value >= (total_marker * lineage_decision_threshold):
                lineage_final_result_dict[key] = value

    lineage_final_result_sorted_dict = OrderedDict(sorted(lineage_final_result_dict.items())) ## sort dict by key
#########################################################

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


