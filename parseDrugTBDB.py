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

    Additional feature HGVS to Genomic positon converter (Beta test) 
    Convert hgvs mutation to genome position based on GFF3 file
    Currently avalible only Substitute,del,ins,dub,inv
    Convert to additional txt file with option -c (indicate convert output filename)
    (Noted: Need deep checking for convertion correctness)
"""

import sys
import json
import argparse
import os
import re

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

amino_code = {"Phe":["UUU","UUC"],
                "Leu":["UUA","UUG","CUU","CUC","CUA","CUG"],
                "Ser":["UCU","UCC","UCS","UCG","AGU","AGC"],
                "Tyr":["UAU","UAC"],
                "Stop":["UAA","UAG","UGA"],
                "Cys":["UGU","UGC"],
                "Trp":["UGG"],
                "Pro":["CCU","CCC","CCA","CCG"],
                "His":["CAU","CAC"],
                "Gln":["CAA","CAG"],
                "Arg":["CGU","CGC","CGA","CGG","AGA","AGG"],
                "Ile":["AUU","AUC","AUA"],
                "Met":["AUG"],
                "Thr":["ACU","ACC","ACA","ACG"],
                "Asn":["AAU","AAC"],
                "Lys":["AAA","AAG"],
                "Val":["GUU","GUC","GUA","GUG"],
                "Ala":["GCU","GCC","GCA","GCG"],
                "Asp":["GAU","GAC"],
                "Glu":["GAA","GAG"],
                "Gly":["GGU","GGC","GGA","GGG"]}

parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i', '--input', help='Input drug db text file.', required=True)
parser.add_argument('-o', '--output', help='Output drug db json file', required=True)
parser.add_argument('-c', '--outputcov', help='Output convert genomic position file', required=True)
parser.add_argument('-g', '--gff3', help='gff3 file for create genomic position file', required=True)

args = parser.parse_args()

input_file = args.input
output_file = args.output
gff3_file = args.gff3
gcov_file = args.outputcov

def readGFF3collectGeneinfo(inputGFF3):
    
    gene_dict = dict() # contain gene and it start stop
    
    with open(inputGFF3,"r",encoding='utf8') as f:
        for line in f:
            no_new_line = line.splitlines()
            info = no_new_line[0].split("\t")
            if info[0][0] == "#":
                continue
            
            chr = info[0]
            db = info[1]
            biotype = info[2]
            start = int(info[3])
            stop = int(info[4])
            strand = info[6]
            bioinfo = info[8]

            if biotype == "gene":
                split_bioinfo = bioinfo.split(";")
                # add gene name info
                gene_name = split_bioinfo[2].split("=")[1]
                if gene_name in gene_dict:
                    print("impossible gene not unique")
                else:
                    gene_dict[gene_name] = [start,stop]

                # add RV ID info
                rv_id = split_bioinfo[0].split("-")[1]
                if rv_id in gene_dict:
                    print("impossible gene not unique")
                else:
                    gene_dict[rv_id] = [start,stop]

    return gene_dict


def resolveGenomePosForProteinHGVS(proteinHGVS,gene_start_pos, protein_pos):

    genome_pos_list = list() ## list contain adjust genome positon. It represent range of two postion. So, finally it must contain only two position. 

    ## extract only sub string from hgvs code
    extract_substring = ",".join(re.findall("[a-zA-Z]+", proteinHGVS))
    extract_substring_list = extract_substring.split(",")

    if "*" in proteinHGVS:
        ## change to stop codon case
        old_amino_acid = extract_substring_list[1]
        new_amino_acid = "Stop"

        old_amino_code_list = amino_code[old_amino_acid]
        new_amino_code_list = amino_code[new_amino_acid]
        all_diff_position_list = list()

        # Loop to collect diff position betweed amino acid code of two amino acid
        one_diff_pos_list = list()
        two_diff_pos_list = list()
        three_diff_pos_list = list()
        count_1_flag = False
        count_2_flag = False
        count_3_flag = False
        for old_amino_code in old_amino_code_list:
            for new_amino_code in new_amino_code_list:
                ## compare two string return diff position (for resolve exact genomic position for p. prefix type)
                ## right now we use all three position (update this later) 
                diff_position_list = [i for i in range(len(old_amino_code)) if old_amino_code[i] != new_amino_code[i]]
                ##############
                if len(diff_position_list) == 1:
                    count_1_flag = True
                    one_diff_pos_list = one_diff_pos_list + diff_position_list
                elif len(diff_position_list) == 2:
                    count_2_flag = True
                    two_diff_pos_list = two_diff_pos_list + diff_position_list
                elif len(diff_position_list) == 3:
                    count_3_flag = True
                    three_diff_pos_list = three_diff_pos_list + diff_position_list

        # select diff pos list. By checking flag in this following order
        if count_1_flag == True:
            all_diff_position_list = one_diff_pos_list
        if count_2_flag == True:
            all_diff_position_list = two_diff_pos_list
        if count_3_flag == True:
            all_diff_position_list = three_diff_pos_list

        # Resolve genome position with adjust position derive from unique diff position
        unique_diff_position_list = list(set(all_diff_position_list))
        for diff_pos in unique_diff_position_list:
            adjust_posiiton = 2 - diff_pos  # reverse adjust number eg diff pos at 0 idx mean adjust number should be 2
            genome_pos = ((gene_start_pos + (protein_pos*3)) - 1) - adjust_posiiton
            genome_pos_list.append(genome_pos)
        
        genome_pos_list.sort()
        #if len(genome_pos_list) > 2:
            ## pick only first and last value. Because final genome_pos_list are represent in range betweet two number         
        #    return [genome_pos_list[0],genome_pos_list[-1]]
        #else:
        genome_pos_string = ",".join(map(str,genome_pos_list))
        return genome_pos_string


    elif "del" not in proteinHGVS and "ins" not in proteinHGVS and "dup" not in proteinHGVS:
        ## Substitution varaint case
        old_amino_acid = extract_substring_list[1]
        new_amino_acid = extract_substring_list[2]

        old_amino_code_list = amino_code[old_amino_acid]
        new_amino_code_list = amino_code[new_amino_acid]

        all_diff_position_list = list()

        # Loop to collect diff position betweed amino acid code of two amino acid
        one_diff_pos_list = list()
        two_diff_pos_list = list()
        three_diff_pos_list = list()
        count_1_flag = False
        count_2_flag = False
        count_3_flag = False
        for old_amino_code in old_amino_code_list:
            for new_amino_code in new_amino_code_list:
                ## compare two string return diff position (for resolve exact genomic position for p. prefix type)
                ## right now we use all three position (update this later) 
                diff_position_list = [i for i in range(len(old_amino_code)) if old_amino_code[i] != new_amino_code[i]]
                ##############
                if len(diff_position_list) == 1:
                    count_1_flag = True
                    one_diff_pos_list = one_diff_pos_list + diff_position_list
                elif len(diff_position_list) == 2:
                    count_2_flag = True
                    two_diff_pos_list = two_diff_pos_list + diff_position_list
                elif len(diff_position_list) == 3:
                    count_3_flag = True
                    three_diff_pos_list = three_diff_pos_list + diff_position_list
        
        # select diff pos list. By checking flag in this following order
        if count_1_flag == True:
            all_diff_position_list = one_diff_pos_list
        if count_2_flag == True:
            all_diff_position_list = two_diff_pos_list
        if count_3_flag == True:
            all_diff_position_list = three_diff_pos_list

        # Resolve genome position with adjust position derive from unique diff position
        unique_diff_position_list = list(set(all_diff_position_list))
        for diff_pos in unique_diff_position_list:
            adjust_posiiton = 2 - diff_pos  # reverse adjust number eg diff pos at 0 idx mean adjust number should be 2
            nucleotide_pos = (protein_pos*3)
            genome_pos = ((gene_start_pos + nucleotide_pos) - 1) - adjust_posiiton
            genome_pos_list.append(genome_pos)
        
        genome_pos_list.sort()
        #if len(genome_pos_list) > 2:
            ## pick only first and last value. Because final genome_pos_list are represent in range betweet two number         
        #    return [genome_pos_list[0],genome_pos_list[-1]]
        #else:
        genome_pos_string = ",".join(map(str,genome_pos_list))
        return genome_pos_string
    
    else:
        ## multiple allele variant case. No resolve position needed
        genome_pos = (gene_start_pos + (protein_pos*3)) - 1
        genome_pos_list.append(genome_pos-2)
        genome_pos_list.append(genome_pos)
        genome_pos_string = str(genome_pos-2)+","+str(genome_pos-1)+","+str(genome_pos)
    
    return genome_pos_string

def convertGenomicPosition(gene_dict,gene_name,hgvs_mutation):
    genome_pos_list = []
    if gene_name in gene_dict:
        gene_start = gene_dict[gene_name][0]
        gene_stop = gene_dict[gene_name][1]

        prefix = hgvs_mutation.split(".")[0]

        ## use re library to extract only prositive and negative number from string
        raw_position = [int(d) for d in re.findall(r'-?\d+', hgvs_mutation)]

        #if len(raw_position) == 0:
            #genome_pos_list.append("0")

        #check prefix and do convertion specific for prefix
        if prefix == "c" or prefix == "r" or prefix == "n" or prefix == "m":
            for pos in raw_position:
                # check positive or negative
                if pos > 0:
                    genome_pos = (gene_start + pos) - 1
                else:
                    genome_pos = (gene_start - pos)
                genome_pos_list.append(genome_pos)
        elif prefix == "p":
            for pos in raw_position:
                # check positive or negative
                if pos > 0:
                    resolve_genome_pos = resolveGenomePosForProteinHGVS(hgvs_mutation,gene_start,pos)
                    genome_pos_list.append(resolve_genome_pos)
                    #genome_pos = (gene_start + (pos*3)) - 1
                    #genome_pos_list.append(genome_pos-2)
                    #genome_pos_list.append(genome_pos)
                else:
                    genome_pos = "protein pos got minus value"
                    genome_pos_list.append(genome_pos)
        elif prefix == "g":
            for pos in raw_position:
                # check positive or negative
                if pos > 0:
                    genome_pos_list.append(pos)
                else:
                    genome_pos = "genomic pos got minus value"
                    genome_pos_list.append(genome_pos)
        else:
            # not hgvs format
            genome_pos_list.append("0")

    else:
        genome_pos_list.append("0")

    result = map(str,genome_pos_list)
    return result


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
gff3_gene_dict = readGFF3collectGeneinfo(gff3_file)
drug_resist_db = dict()
header = False

gcov = open(gcov_file, "w")

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
            not_protein_flag = 0
            # Recorrect mutation format [Recorrect for WHO mutatuon symbol]
            if len(mutation.split(".")) == 1 and len(mutation.split(" ")) == 1 and mutation != "frameshift":
                dummy = ""
                for mutation_str in mutation:
                    if mutation_str == "(":
                        continue
                    elif mutation_str == ")":
                        continue
                    elif mutation_str in protein_code:
                        dummy = dummy + protein_code[mutation_str]
                        not_protein_flag = 1
                    else:
                        dummy = dummy + mutation_str
                if not_protein_flag == 1: ## if it not protien. it is special case. we will not change it format.
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

            #######
            # convert hgvs mutation to genome positon
            # create separate file indicate genoome posiiton and hgvs
            convert_res = convertGenomicPosition(gff3_gene_dict,gene,mutation)
            gcov.write(gene + "\t" + mutation + "\t" + '-'.join(convert_res))
            gcov.write("\n")
            ##########################################################

        header = False

gcov.close()

with open(output_file, 'w') as fp:
    json.dump(drug_resist_db, fp)









