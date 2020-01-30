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

input_folder = ""


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
    drug_typ_itol_path = os.path.join(out_dir, "drug_typ.itol.txt")

    lineage_itol = open(lineage_itol_path, 'w')
    drug_itol = open(drug_itol_path, 'w')
    drug_typ_itol = open(drug_typ_itol_path, 'w')

    #### initiate header of lineage itol anotation file ####

    lineage_itol.write("DATASET_COLORSTRIP\n")
    lineage_itol.write("SEPARATOR\tTAB\n")
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
    drug_itol.write("SEPARATOR\tTAB\n")
    drug_itol.write("DATASET_LABEL\tDrugs\n")
    drug_itol.write("COLOR\t#ff0000\n\n")

    drug_itol.write("SHOW_LEBELS\t1\n")
    drug_itol.write("FIELD_SHAPES\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\t2\n")
    drug_itol.write("FIELD_COLORS\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\tblack\n")
    drug_itol.write("FIELD_LABELS\trifampicin\tisoniazid\tethambutol\tpyrazinamide\tstreptomycin\tfluoroquinolones\taminoglycosides\tkanamycin\tamikacin\tcapreomycin\tethionamide\tpara-aminosalicylic_acid\tclofazimine\tlinezolid\tbedaquiline\n\n")

    drug_itol.write("DATA\n")

    drug_template_map = {"rifampicin":"0","isoniazid":"0","ethambutol":"0","pyrazinamide":"0","streptomycin":"0","fluoroquinolones":"0","aminoglycosides":"0","kanamycin":"0","amikacin":"0","capreomycin":"0","ethionamide":"0","para-aminosalicylic_acid":"0","clofazimine":"0","linezolid":"0","bedaquiline":"0"}
    ################################################

    #### initiate header of drug resist type itol anotation file ####

    drug_typ_itol.write("DATASET_COLORSTRIP\n")
    drug_typ_itol.write("SEPARATOR\tTAB\n")
    drug_typ_itol.write("DATASET_LABEL\tDrug-Resistant\n")
    drug_typ_itol.write("COLOR\t#ff0000\n\n")

    drug_typ_itol.write("LEGEND_TITLE\tDrug resistance\n")
    drug_typ_itol.write("LEGEND_SHAPES\t1\t1\t1\t1\n")
    drug_typ_itol.write("LEGEND_COLORS\t#80FF00\t#7fe5f0\t#8000FF\t#FF0000\n")
    drug_typ_itol.write("LEGEND_LABELS\tSensitive\tDrug-resisant\tMDR\tXDR\n\n")

    drug_typ_itol.write("DATA\n")

    drug_typ_color_map = {"sensitive":"#80FF00","resistant":"#7fe5f0","MDR":"#8000FF","XDR":"#FF0000"}
    ################################################

    list_file = get_in_file_list(in_dir, file_target)
    for file in list_file:
        with open(file) as fp:
            data_dict = json.load(fp)

            #### Extract Data ####

            master_lineage_list = dict()
            lineage_dict = data_dict["lineage"]

            for lineage in lineage_dict:
                #### convert lineageBOV and lineageBOV_AFRI to bovis then check main first then check for mix infect


            master_drug_dict = dict()
            drug_resist_type = data_dict["drug_resist_type"]
            drug_resist_dict = data_dict["small_variant_dr"]
            drug_resist_list = list()
            for key in drug_resist_dict:
                info_variant_dict = drug_resist_dict[key]
                drug_info_dict = info_variant_dict["Drug"]

                for drug in drug_info_dict:
                    if drug not in master_drug_dict:
                        confident = drug_info_dict[drug]
                        master_drug_dict[drug] = confident
            #######################

            #### create itol lineage annotation file ####



            #############################################


            lineage_list = data_dict["lineage"]



f_list = get_in_file_list(in_dir, input_file_extension)