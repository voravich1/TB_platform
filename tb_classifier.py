import json
import vcfpy
from pathlib import Path

vcfFile_name = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/97_typing_snp_vep_custom_ann_hgvs_vcftest.vcf"
drugDBFile_name = "/Users/worawich/Download_dataset/TB_platform_test/drug_db/tbprofiler_drugDB.json"

reader = vcfpy.Reader.from_path(vcfFile_name)

#vcfpy.record.Record.INFO

header = reader.header
vep_header = header.get_info_field_info('CSQ')
vep_header_description = vep_header.description
vep_field_order = vep_header_description.split(': ')[1].split('|')
vep_field_query_index = dict()
hgvsc_idx = 0
hgvsp_idx = 0
count = 0

hgvsc_dict = dict()
hgvsp_dict = dict()

drugDB_gene_dict = dict()
drugDB_hgvs_dict = dict()
drugDB_drug_dict = dict()

result_dict = dict()

###########
## Loop check hgvs field index
###########

for vep_field in vep_field_order:
    vep_field_query_index[vep_field] = count

    if(vep_field == 'HGVSc'):
        hgvsc_idx = count
    elif(vep_field == 'HGVSp'):
        hgvsp_idx = count

    count+=1


##########
# Loop file drug database
# drug database of TBprofiler is json format
# transform DB structure to NBT_AON structure [Gene => HGVS => Drug => Confident]
##########

with open(drugDBFile_name) as json_file:
    data = json.load(json_file)

    for gene in data:

        gene_dict = data[gene]

        drugDB_hgvs_dict = dict()
        for variant in gene_dict:

            variant_dict = gene_dict[variant]

            drug_dict = variant_dict['drugs']
            hgvs = variant_dict['hgvs_mutation']

            drugDB_drug_dict = dict()
            for drug_name in drug_dict:

                drug_name_dict = drug_dict[drug_name]
                confidence = drug_name_dict['confidence']

                drugDB_drug_dict[drug_name] = confidence

            drugDB_hgvs_dict[hgvs] = drugDB_drug_dict

        drugDB_gene_dict[gene] = drugDB_hgvs_dict


###########
## Loop Filter and Collect only record that give HGVS field
## (My understanding right now is "record that has HGVSc and HGVSp is variant record that occur on gene")
###########

record_hit_ID = 0

for record in reader:

    info = record.INFO
    record_csq = info.get('CSQ')


    for vep_field in record_csq:
        vep_field_list = vep_field.split("|")

        if(vep_field_list[hgvsc_idx] == "" and vep_field_list[hgvsp_idx] == ""):
            continue

        vep_gene = vep_field_list[vep_field_query_index['Gene']]

        hgvsc_info = vep_field_list[hgvsc_idx]
        hgvsp_info = vep_field_list[hgvsp_idx]

        if hgvsp_info != "":
            hgvsc_ready = hgvsc_info.split(':', 1)[1]
        else:
            hgvsc_ready = hgvsc_info

        if hgvsp_info != "":
            hgvsp_ready = hgvsp_info.split(':', 1)[1]
        else:
            hgvsp_ready = hgvsp_info

        # check with database
        if vep_gene in drugDB_gene_dict:
            gene_hit_dict = drugDB_gene_dict[vep_gene]

            record_hit_dict = dict()
            if hgvsc_ready in gene_hit_dict:

                record_hit_dict['Drug'] = gene_hit_dict
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(vep_field_order)):
                    record_hit_dict[vep_field_order[i]] = vep_field_list[i]

                result_dict[record_hit_ID] = record_hit_dict

            elif hgvsp_ready in gene_hit_dict:

                record_hit_dict['Drug'] = gene_hit_dict
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(vep_field_order)):
                    record_hit_dict[vep_field_order[i]] = vep_field_list[i]

                result_dict[record_hit_ID] = record_hit_dict

        if hgvsc_ready in hgvsc_dict:
            print("repeat!!" + hgvsc_ready)
            record_list = hgvsc_dict[hgvsc_ready]
            record_list.append(record)
            hgvsc_dict[hgvsc_ready] = record_list
        else:
            record_list = [record]
            hgvsc_dict[hgvsc_ready] = record_list

        if hgvsp_ready in hgvsp_dict:
            print("repeat!!" + hgvsp_ready)
            record_list = hgvsp_dict[hgvsp_ready]
            record_list.append(record)
            hgvsp_dict[hgvsp_ready] = record_list
        else:
            record_list = [record]
            hgvsp_dict[hgvsp_ready] = record_list

print("Done read vcf")






