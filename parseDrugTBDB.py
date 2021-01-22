#!/usr/bin/env python

"""This program is use to parse drug database from tab delimit (txt) format to json format
    Compatible for use with tbprofiler and aonprofiler

    Capable to convert single character protein variant annotation to
    three character protein varaint annotation (HGVS format)
    It best to have mutation in HGVS format

    Did not do anything on other annotation format
    For other format please manually change it by your self. Sorry +_+"

    Input txt format must contain 5 column
    1. Drug - Name of the drug.
    2. Gene - Can be gene names (e.g. inhA, katG) or locus tag (e.g. Rv0678, Rv0682)
    3. Mutation - These must be hgvs nomenclature (e.g. p.Val139Leu, c.-15C>T)
    4. Resistance Level - A resistance level base on MIC level or other testing. Can specify into four level “unknown, Low, Moderate and High”
    5. Confidence Level - These contain text indicate how confidence of this mutation to be drug resistance mutation. Can specify into four level “unknown, Low, Moderate and High”
    Noted: All gene and mutation must based on H37RV reference

"""

import sys
import json
import argparse
import os

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


parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i', '--input', help='Input drug db text file.', required=True)
parser.add_argument('-o', '--output', help='Output drug db json file', required=True)

args = parser.parse_args()

input_file = args.input
output_file = args.output


def compare_level(old_res_level,new_res_level,old_con_level,new_con_level):

    res_level_dict = {"high":3,"moderate":2,"low":1,"unknown":0}
    con_level_dict = {"high":3,"moderate":2,"low":1,"indeterminate":0}

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


#tf = open(tranform_input_file,"w")

drug_resist_db = dict()
header = False

with open(input_file,"r",encoding='utf8') as f:
    for line in f:
        no_new_line = line.splitlines()
        info = no_new_line[0].split("\t")
        drug = info[0].lower()
        gene = info[1]
        mutation = info[2]
        resist_level = info[3].lower()
        confidence_level = info[4].lower()

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

            if dummy_resist_level == "high":
                resist_level = "high"
            elif dummy_resist_level == "moderate":
                resist_level = "moderate"
            elif dummy_resist_level == "low":
                resist_level = "low"
            else:
                resist_level = "unknown"
            #################################

            # Recorrect confidence
            if confidence_level == "high":
                confidence_level = "high"
            elif confidence_level == "moderate":
                confidence_level = "moderate"
            elif confidence_level == "low":
                confidence_level = "low"
            else:
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

        header = False


with open(output_file, 'w') as fp:
    json.dump(drug_resist_db, fp)









