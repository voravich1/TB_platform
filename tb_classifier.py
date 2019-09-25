import json
import vcfpy
from pathlib import Path

vcfFile_name = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/97_typing_snp_vep_custom_ann_hgvs_vcftest.vcf"

reader = vcfpy.Reader.from_path(vcfFile_name)

#vcfpy.record.Record.INFO

header = reader.header
vep_header = header.get_info_field_info('CSQ')
vep_header_description = vep_header.description
vep_field_order = vep_header_description.split(': ')[1].split('|')
hgvsc_idx = 0
hgvsp_idx = 0
count = 0

###########
## Loop check hgvs field index
###########

for vep_field in vep_field_order:
    if(vep_field == 'HGVSc'):
        hgvsc_idx = count
    elif(vep_field == 'HGVSp'):
        hgvsp_idx = count

    count+=1
    vcfpy.Record
###########
## Loop Filter and Collect only record that give HGVS field
## (My understanding right now is "record that has HGVSc and HGVSp is variant record that occur on gene")
###########

for record in reader:

    info = record.INFO
    csq = info.get('CSQ')


    for vep_field in csq:
        vep_field_list = vep_field.split("|")

        if(vep_field_list[hgvsc_idx] == "" and vep_field_list[hgvsp_idx] == ""):
            continue



    print("")
