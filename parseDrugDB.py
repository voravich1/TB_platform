#!/usr/bin/env python

"""This program is use to parse drug database from normal format to json format
    Compatible for use with tbprofiler and aonprofiler

    Focus on convert single character protein variant annotation to
    three character protein varaint annotation (HGVS format)

    Did not do anything on other annotation format
    For other format please manually change it by your self
"""

import sys
import json
import csv
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

protein_code = {"A":"Ala","R":"Arg","N":"Asn","D":"Asp","B":"Asx",
                "C":"Cys","E":"Glu","Q":"Gln","Z":"Glx","G":"Gly",
                "H":"His","I":"Ile","L":"Leu","K":"Lys","M":"Met",
                "F":"Phe","P":"Pro","S":"Ser","T":"Thr","W":"Trp",
                "Y":"Tyr","V":"Val"}

input_file = "/Users/worawich/Download_dataset/tb_platform/drug_db_committee/thai_drug_db_common.csv" # csv file contain drug resistant information
output_file = "/Users/worawich/Download_dataset/tb_platform/drug_db_committee/thai_drug_db.json"
tranform_input_file = "/Users/worawich/Download_dataset/tb_platform/drug_db_committee/TB_thai_drug_db.csv"

def compare_level(old_res_level,new_res_level,old_con_level,new_con_level):

    res_level_dict = {"high":3,"moderate":2,"low":1,"unknown":0}
    con_level_dict = {"high":3,"moderate":2,"minimal":1,"indeterminate":0}

    ## resistant level compare
    if res_level_dict[old_res_level] > res_level_dict[new_res_level]:
        select_res_level =  old_res_level
    elif res_level_dict[old_res_level] < res_level_dict[new_res_level]:
        select_res_level = new_res_level
    else:
        select_res_level = old_res_level

    ## confidence level compare
    if con_level_dict[old_con_level] > con_level_dict[new_con_level]:
        select_con_level = old_con_level
    elif con_level_dict[old_con_level] < con_level_dict[new_con_level]:
        select_con_level = new_con_level
    else:
        select_con_level = old_con_level

    return select_res_level, select_con_level


tf = open(tranform_input_file,"w")

drug_resist_db = dict()
header = True

with open(input_file,"r",encoding='utf8') as f:
    for line in f:
        no_new_line = line.splitlines()
        info = no_new_line[0].split(",")
        drug = info[0]
        gene = info[1]
        mutation = info[2]
        resist_level = info[3]
        confidence_level = info[4]

        if header != True :

            # Recorrect mutation format
            if len(mutation.split(".")) == 1 and len(mutation.split(" ")) == 1 and mutation != "frameshift":
                dummy = ""
                for str in mutation:
                    if str == "(":
                        continue
                    elif str == ")":
                        continue
                    elif str in protein_code:
                        dummy = dummy + protein_code[str]
                    else:
                        dummy = dummy + str
                mutation = "p." + dummy
            #################################

            # Recorrect resistant level
            if len(resist_level.split("(")) > 1:
                dummy_resist_level = resist_level.split("(")[0]
            else:
                dummy_resist_level = resist_level

            if dummy_resist_level == "High MIC mutation":
                resist_level = "high"
            elif dummy_resist_level == "Moderate MIC mutation":
                resist_level = "moderate"
            elif dummy_resist_level == "Low MIC mutation":
                resist_level = "low"
            elif dummy_resist_level == "Resistance level unknown":
                resist_level = "unknown"
            #################################

            # Recorrect confidence
            if confidence_level == "High":
                confidence_level = "high"
            elif confidence_level == "Moderate":
                confidence_level = "moderate"
            elif confidence_level == "Minimal":
                confidence_level = "minimal"
            elif confidence_level == "-":
                confidence_level = "indeterminate"
            #################################


            # put data to drug db dict
            if gene in drug_resist_db:
                variant_dict = drug_resist_db[gene]

                if mutation in variant_dict:
                    inner_variant_dict = variant_dict[mutation]
                    drug_dict = inner_variant_dict["drugs"]
                    inner_variant_dict["hgvs_mutation"] = mutation
                    if drug in drug_dict:
                        print("impossible\n")
                        confidence_resist_level_dict = drug_dict[drug]
                        dummy_res_level = confidence_resist_level_dict["resistance level"]
                        dummy_con_level = confidence_resist_level_dict["confidence"]

                        select_res_level, select_con_level = compare_level(dummy_res_level,resist_level,dummy_con_level,confidence_level)

                        confidence_resist_level_dict = {"resistance level": select_res_level, "confidence": select_con_level}
                        drug_dict[drug] = confidence_resist_level_dict

                    else:
                        confidence_resist_level_dict = {"resistance level": resist_level, "confidence": confidence_level}
                        drug_dict[drug] = confidence_resist_level_dict

                    inner_variant_dict["drugs"] = drug_dict
                    variant_dict[mutation] = inner_variant_dict
                else:
                    confidence_resist_level_dict = {"resistance level": resist_level, "confidence": confidence_level}
                    drug_dict = {drug: confidence_resist_level_dict}
                    inner_variant_dict = {"drugs":drug_dict,"hgvs_mutation":mutation}
                    variant_dict[mutation] = inner_variant_dict

                drug_resist_db[gene] = variant_dict
            else:
                confidence_resist_level_dict = {"resistance level":resist_level,"confidence":confidence_level}
                drug_dict = {drug:confidence_resist_level_dict}
                inner_variant_dict = {"drugs":drug_dict,"hgvs_mutation":mutation}
                variant_dict = {mutation:inner_variant_dict}
                drug_resist_db[gene] = variant_dict
        tf.write(drug + "," + gene + "," + mutation + "," + resist_level + "," + confidence_level + "\n")
        header = False
tf.close()

with open(output_file, 'w') as fp:
    json.dump(drug_resist_db, fp)









