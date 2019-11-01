import json
import vcfpy
from collections import OrderedDict
from pathlib import Path

vcfFile_vep_name = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/97_typing_snp_vep_custom_ann_hgvs_vcftest.vcf"

vcfFile_snpeff_name = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/97_typing_snp_snpeff.vcf"

drugDBFile_name = "/Users/worawich/Download_dataset/TB_platform_test/drug_db/tbprofiler_drugDB.json"

json_result_file = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/97_lineage_drug_result2.json"

reader_vep_vcf = vcfpy.Reader.from_path(vcfFile_vep_name)
reader_snpeff_vcf = vcfpy.Reader.from_path(vcfFile_snpeff_name)

header_vep_vcf = reader_vep_vcf.header
vep_header = header_vep_vcf.get_info_field_info('CSQ')
vep_header_description = vep_header.description
vep_field_order = vep_header_description.split(': ')[1].split('|')
vep_field_query_index = dict()
hgvsc_vep_idx = 0
hgvsp_vep_idx = 0
hgvsc_vep_dict = dict()
hgvsp_vep_dict = dict()

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
## Loop check hgvs field index in vep vcf
###########
count = 0
for vep_field in vep_field_order:
    vep_field_query_index[vep_field] = count
    if(vep_field == 'HGVSc'):
        hgvsc_vep_idx = count
    elif(vep_field == 'HGVSp'):
        hgvsp_vep_idx = count

    count+=1

#################################################

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

            if variant == '-15C>T':
                print("")

            variant_dict = gene_dict[variant]

            drug_dict = variant_dict['drugs']
            hgvs = variant_dict['hgvs_mutation']

            if hgvs == "p.Leu452Pro":
                print()

            drugDB_drug_dict = dict()
            for drug_name in drug_dict:

                drug_name_dict = drug_dict[drug_name]
                confidence = drug_name_dict['confidence']

                drugDB_drug_dict[drug_name] = confidence

            drugDB_hgvs_dict[hgvs] = drugDB_drug_dict

        drugDB_gene_dict[gene] = drugDB_hgvs_dict

###########################################################


###########
## For snpEff_VCF
## Loop Filter and Collect only record that give HGVS field
## Goal for loop on snpEff VCF is to get drug information
## (My understanding right now is "record that has HGVSc and HGVSp is variant record that occur on gene")
###########

record_hit_ID = 0

for record in reader_snpeff_vcf:

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

        if hgvsp_ready == "p.Leu452Pro":
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

print("Done read vcf")

###########
## For VEP_VCF
## Loop Filter and Collect only record that give HGVS field
## Goal for loop on VEP VCF is to get lineage information
## (My understanding right now is "record that has HGVSc and HGVSp is variant record that occur on gene")
###########

record_hit_ID = 0

for record in reader_vep_vcf:

    info = record.INFO
    record_csq = info.get('CSQ')


    for vep_field in record_csq:
        vep_field_list = vep_field.split("|")

        if(vep_field_list[hgvsc_vep_idx] == "" and vep_field_list[hgvsp_vep_idx] == ""):
            continue

        vep_gene = vep_field_list[vep_field_query_index['Gene']]
        #vep_lineage = vep_field_list[vep_field_query_index['TB_DB_LIN']]
        vep_lineage = vep_field_list[len(vep_field_list)-1]

        ####################################################
        ## Count lineage
        ####################################################
        if vep_lineage != "":
            if vep_lineage in lineage_result_dict:
                lineage_count = lineage_result_dict[vep_lineage]
                lineage_count += 1
                lineage_result_dict[vep_lineage] = lineage_count
            else:
                lineage_result_dict[vep_lineage] = 1

        ######################################################

        hgvsc_info = vep_field_list[hgvsc_vep_idx]
        hgvsp_info = vep_field_list[hgvsp_vep_idx]

        if hgvsp_info != "":
            hgvsc_ready = hgvsc_info.split(':', 1)[1]
        else:
            hgvsc_ready = hgvsc_info

        if hgvsp_info != "":
            hgvsp_ready = hgvsp_info.split(':', 1)[1]
        else:
            hgvsp_ready = hgvsp_info

       # if hgvsp_ready == "p.Leu452Pro":
        #    print("")

        #if hgvsc_ready == 'c.-15C>T':
         #   print(vep_gene)

        #if hgvsp_ready == 'c.-15C>T':
         #   print(vep_gene)

        # check with database
        if vep_gene in drugDB_gene_dict:
            gene_hit_dict = drugDB_gene_dict[vep_gene]

            record_hit_dict = dict()
            if hgvsc_ready in gene_hit_dict:

                record_hit_dict['Drug'] = gene_hit_dict[hgvsc_ready]
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(vep_field_order)):
                    record_hit_dict[vep_field_order[i]] = vep_field_list[i]

                drug_result_dict[record_hit_ID] = record_hit_dict
                record_hit_ID += 1

            elif hgvsp_ready in gene_hit_dict:

                record_hit_dict['Drug'] = gene_hit_dict[hgvsp_ready]
                record_hit_dict['Ref'] = record.REF
                record_hit_dict['Pos'] = record.POS

                for i in range(len(vep_field_order)):
                    record_hit_dict[vep_field_order[i]] = vep_field_list[i]

                drug_result_dict[record_hit_ID] = record_hit_dict
                record_hit_ID += 1

        if hgvsc_ready in hgvsc_vep_dict:
            #print("repeat!!" + hgvsc_ready)
            record_list = hgvsc_vep_dict[hgvsc_ready]
            record_list.append(record)
            hgvsc_vep_dict[hgvsc_ready] = record_list
        else:
            record_list = [record]
            hgvsc_vep_dict[hgvsc_ready] = record_list

        if hgvsp_ready in hgvsp_vep_dict:
            #print("repeat!!" + hgvsp_ready)
            record_list = hgvsp_vep_dict[hgvsp_ready]
            record_list.append(record)
            hgvsp_vep_dict[hgvsp_ready] = record_list
        else:
            record_list = [record]
            hgvsp_vep_dict[hgvsp_ready] = record_list

#########################################################


######################################
## Select majority vote on lineage result
######################################
itemMaxValue = max(lineage_result_dict.items(), key=lambda x: x[1])


listOfMaxKeys = list()
# Iterate over all the items in dictionary to find keys with max value
for key, value in lineage_result_dict.items():
    if value == itemMaxValue[1]:
        listOfMaxKeys.append(key)

#print('Keys with maximum Value in Dictionary : ', listOfKeys)

for lineage in lineage_result_dict:
    main_lineage = lineage.split(".")[0]
    if main_lineage in listOfMaxKeys:
        lineage_final_result_dict[lineage] = lineage_result_dict[lineage]


lineage_final_result_sorted_dict = OrderedDict(sorted(lineage_final_result_dict.items()))


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


if ('isoniazid' in drug_all_hit) and ("rifampicin" in drug_all_hit):

    # check fluoroquinolone derivative
    if ('ciprofloxacin' in drug_all_hit) or ('Garenoxacin'  in drug_all_hit) or ('Gatifloxacin' in drug_all_hit) or ('Gemifloxacin' in drug_all_hit) or ('Levofloxacin' in drug_all_hit) or ('Moxifloxacin' in drug_all_hit):

        # check with secondline injectable drug
        if ('capreomycin' in drug_all_hit) or ('kanamycin' in drug_all_hit) or ('amikacin' in drug_all_hit):

            drug_resist_type = "XDR"
        else:
            drug_resist_type = "MDR"
    else:
        drug_resist_type = "MDR"
else:
    drug_resist_type = ""
###########################################

###########################################
## Create json dict contain lineage, drug result and drug resistant type
## Then save to json file
###########################################

result_dict["lineage"] = lineage_final_result_sorted_dict
result_dict["small_variant_dr"] = drug_result_dict
result_dict["drug_resist_type"] = drug_resist_type
js = json.dumps(result_dict, sort_keys=True, indent=4)


json_write = open(json_result_file,'w')

json_write.write(js)

json_write.close()
print("Done read vcf")






