#! /Users/worawich/miniconda3/envs/tbprofiler/bin/python

import sys
from typing import List, Any, Dict

import acmg_grader
import json

"""
    This program will grade drug resistance follow the acmg grading criteria
    first input is acmg grading database file
    second input is tbprofiler result txt format
"""
acmg_grading_db = sys.argv[1]
tbprofiler_json = sys.argv[2]
save_json_grade_result = sys.argv[3]

gene_hgvs_map: Dict[str, list] = {}

grader = acmg_grader.Grader(acmg_grading_db)

with open(tbprofiler_json, "r") as data:
    tbprofiler_result = json.load(data)

    drug_variant = tbprofiler_result.get('dr_variants')

    # loop first round for collect total gene hgvs pair
    for variant in drug_variant:
        gene_name = variant.get('gene_name')
        hgvs = variant.get('change')
        drug = variant.get('drug')

        if gene_name in gene_hgvs_map:
            hgvs_list = []
            hgvs_list = gene_hgvs_map.get(gene_name)
            hgvs_list.append(hgvs)
            addData = {gene_name : hgvs_list}
            gene_hgvs_map.update(addData)
        else:
            hgvs_list = []
            hgvs_list.append(hgvs)
            gene_hgvs_map[gene_name] = hgvs_list

    # second loop for grading
    for variant in drug_variant:
        gene_name = variant.get('gene_name')
        hgvs = variant.get('change')
        drug = variant.get('drug')
        grade = grader.getGradeSpecialCase(drug, gene_name, hgvs, gene_hgvs_map)
        addData = {"grade" : grade}
        variant.update(addData)

with open(save_json_grade_result, 'w') as f:
    json.dump(tbprofiler_result, f)




    # check special case


def convertOne2ThreeAmino(oneLetter):

    if(oneLetter == "A"):
        threeLetter="Ala"
        return threeLetter
    elif(oneLetter == "R"):
        threeLetter = "Arg"
        return threeLetter
    elif (oneLetter == "N"):
        threeLetter = "Asn"
        return threeLetter
    elif (oneLetter == "D"):
        threeLetter = "Asp"
        return threeLetter
    elif (oneLetter == "C"):
        threeLetter = "Cys"
        return threeLetter
    elif (oneLetter == "E"):
        threeLetter = "Glu"
        return threeLetter
    elif (oneLetter == "Q"):
        threeLetter = "Gln"
        return threeLetter
    elif (oneLetter == "G"):
        threeLetter = "Gly"
        return threeLetter
    elif (oneLetter == "H"):
        threeLetter = "His"
        return threeLetter
    elif (oneLetter == "I"):
        threeLetter = "Ile"
        return threeLetter
    elif (oneLetter == "L"):
        threeLetter = "Leu"
        return threeLetter
    elif (oneLetter == "K"):
        threeLetter = "Lys"
        return threeLetter
    elif (oneLetter == "M"):
        threeLetter = "Met"
        return threeLetter
    elif (oneLetter == "F"):
        threeLetter = "Phe"
        return threeLetter
    elif (oneLetter == "P"):
        threeLetter = "Pro"
        return threeLetter
    elif (oneLetter == "U"):
        threeLetter = "Gip"
        return threeLetter
    elif (oneLetter == "S"):
        threeLetter = "Ser"
        return threeLetter
    elif (oneLetter == "T"):
        threeLetter = "Thr"
        return threeLetter
    elif (oneLetter == "W"):
        threeLetter = "Trp"
        return threeLetter
    elif (oneLetter == "Y"):
        threeLetter = "Tyr"
        return threeLetter
    elif (oneLetter == "V"):
        threeLetter = "Val"
        return threeLetter




