#!/usr/bin/env python

"""This program is use for collate all json file result from tb_aon_profiler
    Compatible with json file from tb_aon_profiler only
"""

import sys
import os
import sys
import argparse
import json
import vcfpy
from collections import OrderedDict
from pathlib import Path

__author__ = "Worawich Phornsiricharoenphant"
__copyright__ = ""
__credits__ = ["Worawich Phornsiricharonphant"]
__license__ = "GPL-3.0"
__version__ = "1.0"
__maintainer__ = "Worawich Phornsiricharoenphant"
__email__ = "worawich.ph@gmail.com"
__status__ = "Development"


parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i', '--input', help='Input folder name', required=True)
parser.add_argument('-o', '--output', help='Absolute path of Output folder', required=True)
parser.add_argument('-t', '--target', help='Extension of target input file Ex. ".vcf", ".vcf.gz"', required=True)

args = parser.parse_args()


def get_in_file_list(in_dir, file_name_part, path_num=1, alternative_directory=None):
    """Basie function.
    This function is for getting a list of files in the target location based on part of a file name.

    Args:
        in_dir (string): a name of a target folder.
        file_name_part (string): a end-part of the file name that you want to get a list on. another word a target signature of file
        path_num (int): option for using alternative input folder part. 1 is default.
        alternative_directory (string): if your input folder is in a different location
                                        from this script, you can input the path with this parameter.

    Return:
        int: The sample number in the list.

    """
    if path_num == 1:
        cwd = os.path.dirname(os.path.realpath(sys.argv[0]))
    else:
        cwd = alternative_directory
    in_path = os.path.join(cwd, in_dir)
    list_files = [os.path.join(in_path, file) for file in os.listdir(in_path) if file.endswith(file_name_part)]
    return list_files;


def create_itol_lineage_drug_annotation_file(in_dir, out_dir, file_target=".json"):
    lineage_itol_path = os.path.join(out_dir, "lineage.itol.txt")
    drug_itol_path = os.path.join(out_dir, "drug_resist.itol.txt")
    drug_conf_itol_path = os.path.join(out_dir, "drug_resist_confident.itol.txt")
    drug_typ_itol_path = os.path.join(out_dir, "drug_typ.itol.txt")

    lineage_itol = open(lineage_itol_path, 'w')
    drug_itol = open(drug_itol_path, 'w')
    drug_conf_itol = open(drug_conf_itol_path, 'w')
    drug_typ_itol = open(drug_typ_itol_path, 'w')

    #### initiate header of lineage itol anotation file ####

    lineage_itol.write("DATASET_COLORSTRIP\n")
    lineage_itol.write("SEPARATOR TAB\n")
    lineage_itol.write("DATASET_LABEL\tLineage\n")
    lineage_itol.write("COLOR\t#ff0000\n\n")

    lineage_itol.write("LEGEND_TITLE\tLineage\n")
    lineage_itol.write("LEGEND_SHAPES\t1\t1\t1\t1\t1\t1\t1\t1\t1\t1\n")
    lineage_itol.write("LEGEND_COLORS\t#104577\t#ab2323\t#18a68c\t#f68e51\t#7cb5d2\t#fde05e\t#bc94b7\t#f8e0c8\t#808080\t#000000\n")
    lineage_itol.write("LEGEND_LABELS\tLineage1\tLineage2\tLineage3\tLineage4\tLineage5\tLineage6\tLineage7\tBovis\tMixed\tOther\n\n")

    lineage_itol.write("DATA\n")

    lineage_color_map = {"lineage1":"#104577","lineage2":"#ab2323","lineage3":"#18a68c","lineage4":"#f68e51","lineage5":"#7cb5d2","lineage6":"#fde05e","lineage7":"#bc64b7","bovis":"#f8e0c8","mixed":"#808080","other":"#000000"}

    ################################################

    #### initiate header of drug resist itol anotation file ####

    drug_itol.write("DATASET_BINARY\n")
    drug_itol.write("SEPARATOR TAB\n")
    drug_itol.write("DATASET_LABEL\tDrugs\n")
    drug_itol.write("COLOR\t#ff0000\n\n")

    drug_itol.write("SHOW_LABELS\t1\n")
    drug_itol.write("FIELD_SHAPES\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\n")
    drug_itol.write("FIELD_COLORS\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\n")
    drug_itol.write("FIELD_LABELS\trifampicin\tisoniazid\tethambutol\tpyrazinamide\tstreptomycin\tfluoroquinolones\taminoglycosides\tkanamycin\tamikacin\tcapreomycin\tethionamide\tpara-aminosalicylic_acid\tclofazimine\tlinezolid\tbedaquiline\tciprofloxacin\tlevofloxacin\tmoxifloxacin\tofloxacin\n\n")

    drug_itol.write("DATA\n")

    ################################################

    #### initiate header of drug resist type itol anotation file ####

    drug_typ_itol.write("DATASET_COLORSTRIP\n")
    drug_typ_itol.write("SEPARATOR TAB\n")
    drug_typ_itol.write("DATASET_LABEL\tDrug-Resistant\n")
    drug_typ_itol.write("COLOR\t#ff0000\n\n")

    drug_typ_itol.write("LEGEND_TITLE\tDrug resistance\n")
    drug_typ_itol.write("LEGEND_SHAPES\t1\t1\t1\t1\n")
    drug_typ_itol.write("LEGEND_COLORS\t#80FF00\t#7fe5f0\t#8000FF\t#FF0000\n")
    drug_typ_itol.write("LEGEND_LABELS\tSensitive\tDrug-resistant\tMDR\tXDR\n\n")

    drug_typ_itol.write("DATA\n")

    drug_typ_color_map = {"sensitive":"#80FF00","resistant":"#7fe5f0","MDR":"#8000FF","XDR":"#FF0000","none":"#000000"}
    ################################################

    #### initiate header of drug resist confident itol anotation file ####

    drug_conf_itol.write("DATASET_EXTERNALSHAPE\n")
    drug_conf_itol.write("SEPARATOR TAB\n")
    drug_conf_itol.write("DATASET_LABEL\tDrugs with confident\n")
    drug_conf_itol.write("COLOR\t#ff0000\n\n")

    drug_conf_itol.write("SHOW_VALUES\t0\n")
    drug_conf_itol.write("FIELD_SHAPES\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\n")
    drug_conf_itol.write("FIELD_COLORS\t#5E5E5E\t#AC8E68\t#FF453A\t#32D74B\t#0A84FF\t#FF9F0A\t#FFD60A\t#BF5AF2\t#64D2FF\t#3F638B\t#6b3e26\t#007D7D\t#4E8F00\t#011893\t#FF7E79\t#38D4D6\t#431479\t#FF7C00\t#FFADD6\n")
    drug_conf_itol.write("FIELD_LABELS\trifampicin\tisoniazid\tethambutol\tpyrazinamide\tstreptomycin\tfluoroquinolones\taminoglycosides\tkanamycin\tamikacin\tcapreomycin\tethionamide\tpara-aminosalicylic_acid\tclofazimine\tlinezolid\tbedaquiline\tciprofloxacin\tlevofloxacin\tmoxifloxacin\tofloxacin\n\n")

    drug_conf_itol.write("LEGEND_TITLE\tDrugs\n")
    drug_conf_itol.write("LEGEND_SHAPES\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\n")
    drug_conf_itol.write("LEGEND_COLORS\t#5E5E5E\t#AC8E68\t#FF453A\t#32D74B\t#0A84FF\t#FF9F0A\t#FFD60A\t#BF5AF2\t#64D2FF\t#3F638B\t#6b3e26\t#007D7D\t#4E8F00\t#011893\t#FF7E79\t#38D4D6\t#431479\t#FF7C00\t#FFADD6\n")
    drug_conf_itol.write("LEGEND_LABELS\trifampicin\tisoniazid\tethambutol\tpyrazinamide\tstreptomycin\tfluoroquinolones\taminoglycosides\tkanamycin\tamikacin\tcapreomycin\tethionamide\tpara-aminosalicylic_acid\tclofazimine\tlinezolid\tbedaquiline\tciprofloxacin\tlevofloxacin\tmoxifloxacin\tofloxacin\n\n")

    drug_conf_itol.write("DATA\n")

    ################################################

    list_file = get_in_file_list(in_dir, file_target)
    for file in list_file:
        with open(file) as fp:

            data_dict = json.load(fp)

            sample_name = data_dict["sample_name"]

            #### Extract lineage Data ####
            lineage_dict = data_dict["lineage"]
            main_lineage_list = list()
            for lineage in lineage_dict:
                if lineage == "lineageBOV" or lineage == "lineageBOV_AFRI":
                    lineage = "bovis"

                if lineage in lineage_color_map:
                    main_lineage_list.append(lineage)

            if len(main_lineage_list) > 1 :
                final_lineage_judgement = "mixed"
            elif len(main_lineage_list) == 1:
                final_lineage_judgement = main_lineage_list[0]
            else:
                final_lineage_judgement = "other"
            #######################

            #### Export lineage Data to lineage itol annotation file ####
            lineage_color = lineage_color_map[final_lineage_judgement]
            lineage_itol.write(sample_name + "\t" + lineage_color + "\n")
            #######################

            #### Extract drug Data ####


            drug_conf_template_map = {"rifampicin": "0", "isoniazid": "0", "ethambutol": "0", "pyrazinamide": "0",
                                      "streptomycin": "0", "fluoroquinolones": "0", "aminoglycosides": "0",
                                      "kanamycin": "0",
                                      "amikacin": "0", "capreomycin": "0", "ethionamide": "0",
                                      "para-aminosalicylic_acid": "0",
                                      "clofazimine": "0", "linezolid": "0", "bedaquiline": "0", "ciprofloxacin": "0", "levofloxacin": "0", "moxifloxacin": "0", "ofloxacin": "0"}
            drug_template_map = {"rifampicin": "0", "isoniazid": "0", "ethambutol": "0", "pyrazinamide": "0",
                                 "streptomycin": "0", "fluoroquinolones": "0", "aminoglycosides": "0", "kanamycin": "0",
                                 "amikacin": "0", "capreomycin": "0", "ethionamide": "0",
                                 "para-aminosalicylic_acid": "0", "clofazimine": "0", "linezolid": "0",
                                 "bedaquiline": "0", "ciprofloxacin": "0", "levofloxacin": "0", "moxifloxacin": "0", "ofloxacin": "0"}

            master_drug_dict = dict()
            drug_resist_dict = data_dict["small_variant_dr"]
            drug_resist_list = list()
            for key in drug_resist_dict:
                info_variant_dict = drug_resist_dict[key]
                drug_info_dict = info_variant_dict["Drug"]

                for drug in drug_info_dict:

                    if drug in drug_template_map:
                        drug_template_map[drug] = "1"

                    if drug not in master_drug_dict:
                        confidence = drug_info_dict[drug]

                        if confidence == "high":
                            drug_conf_template_map[drug] = "400"
                        elif confidence == "moderate":
                            drug_conf_template_map[drug] = "300"
                        elif confidence == "low":
                            drug_conf_template_map[drug] = "200"
                        elif confidence == "indeterminate":
                            drug_conf_template_map[drug] = "100"
            #######################
            #### Export drug Data to drug itol annotation file ####
            info_list = [sample_name]
            for drug in drug_template_map:
                info_list.append(drug_template_map[drug])
            write_content = '\t'.join(info_list)
            drug_itol.write(write_content + "\n")
            #######################
            #### Export drug with cofidence Data to drug confidence itol annotation file ####
            info_list = [sample_name]
            for drug in drug_conf_template_map:
                info_list.append(drug_conf_template_map[drug])
            write_content = '\t'.join(info_list)
            drug_conf_itol.write(write_content + "\n")
            #######################

            #### Extract drug type Data and Export to drug type itol annotation file ####
            drug_resist_type = data_dict["drug_resist_type"]
            drug_typ_itol.write(sample_name+"\t"+drug_typ_color_map[drug_resist_type]+"\n")
            #######################

    drug_itol.close()
    drug_conf_itol.close()
    drug_typ_itol.close()
    lineage_itol.close()
########################################################################################################################
def collate_result(inputs, out_dir, file_target=args.target):

    collate_file_path = os.path.join(out_dir, "summary_result.txt")
    collate_file = open(collate_file_path, 'w')

    # Initialize header
    collate_file.write("Sample\tLineage\tDrug Resistant Type\trifampicin\tisoniazid\tethambutol\tpyrazinamide\tstreptomycin\tfluoroquinolones\taminoglycosides\tkanamycin\tamikacin\tcapreomycin\tethionamide\tpara-aminosalicylic_acid\tclofazimine\tlinezolid\tbedaquiline\tciprofloxacin\tlevofloxacin\tmoxifloxacin\tofloxacin\n")


    list_file = get_in_file_list(inputs, file_target)
    for file in list_file:
        with open(file) as fp:

            data_dict = json.load(fp)

            sample_name = data_dict["sample_name"]

            # extract lineage data
            lineage_dict = data_dict["lineage"]
            lineage_list = []
            for lineage in lineage_dict:
                lineage_list.append(lineage)
            lineage_export_string = "|".join(lineage_list)

            # extract drug resistant
            drug_template_map = {"rifampicin": 0, "isoniazid": 0, "ethambutol": 0, "pyrazinamide": 0,
                                 "streptomycin": 0, "fluoroquinolones": 0, "aminoglycosides": 0, "kanamycin": 0,
                                 "amikacin": 0, "capreomycin": 0, "ethionamide": 0,
                                 "para-aminosalicylic_acid": 0, "clofazimine": 0, "linezolid": 0,
                                 "bedaquiline": 0, "ciprofloxacin": 0, "levofloxacin": 0, "moxifloxacin": 0,
                                 "ofloxacin": 0}

            drug_resist_dict = data_dict["small_variant_dr"]
            for key in drug_resist_dict:
                info_variant_dict = drug_resist_dict[key]
                drug_info_dict = info_variant_dict["Drug"]

                for drug in drug_info_dict:

                    if drug in drug_template_map:
                        confidence = drug_info_dict[drug]
                        if confidence == "high":
                            score = 4
                        elif confidence == "moderate":
                            score = 3
                        elif confidence == "low":
                            score = 2
                        elif confidence == "indeterminate":
                            score = 1
                        else:
                            score = 0

                        # check existing score value if it lower tand current change to current score
                        if drug_template_map[drug] <= score:
                            drug_template_map[drug] = score
            drug_resist_list = []
            for drug in drug_template_map:
                score = drug_template_map[drug]
                if score == 4:
                    drug_resist_list.append("high")
                elif score == 3:
                    drug_resist_list.append("moderate")
                elif score == 2:
                    drug_resist_list.append("low")
                elif score == 1:
                    drug_resist_list.append("indeterminate")
                else:
                    drug_resist_list.append("none")

            drug_resist_export_string = "\t".join(drug_resist_list)

            # extract drug resistant type
            drug_resist_type = data_dict["drug_resist_type"]
            if drug_resist_type == "":
                drug_resist_type = "none"


            collate_file.write(sample_name + "\t" + lineage_export_string + "\t" + drug_resist_type + "\t" + drug_resist_export_string + "\n")

    collate_file.close()
# --------main---------

#f_list = get_in_file_list(args.input, args.target)
create_itol_lineage_drug_annotation_file(args.input, args.output, file_target=args.target)
collate_result(args.input, args.output, file_target=args.target)
