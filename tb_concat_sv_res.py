import sys
import os
import gzip

######
#
# script for concatinate multiple SV vcf file
# Can be able to use with vcf from delly and manta
#
######

#path="/Volumes/10TBSeagateBackupPlus/NBT/TB_platform/analysis/tb_delly/tb_1188/DellySVResults"
#path="/Volumes/10TBSeagateBackupPlus/NBT/thalassemia/analysis/HS08/analysis/res_delly"
#path="/Volumes/10TBSeagateBackupPlus/NBT/thalassemia/analysis/HS08/analysis/MantaSVvcf"
path="/Volumes/10TBSeagateBackupPlus/NBT/scratchHBA/scratchHBA_region_mantaV1_4_vcf"
#path="/Volumes/4TB_WD/TB/Martin_new_list/TB_result_BWA_SV_Ver/manta/vcf_ann"
list_target_sample_file="/Volumes/4TB_WD/TB/Martin_new_list/TB_result_BWA_SV_Ver/manta/list_L1_1_1.txt"
sum_by_target = False
target_suffix = ".vcf.gz"

list_of_files = {}
if sum_by_target == True :

    target_sample_dict = dict()
    with open(list_target_sample_file, 'r') as list_file:
        for line in list_file:
            samplename = line.splitlines()[0]
            target_sample_dict[samplename] = True

    for (dirpath, dirnames, filenames) in os.walk(path):
        for filename in filenames:
            if filename.endswith(target_suffix):
                if filename in target_sample_dict:
                    list_of_files[filename] = os.sep.join([dirpath, filename])
elif sum_by_target == False:
    for (dirpath, dirnames, filenames) in os.walk(path):
        for filename in filenames:
            if filename.endswith(target_suffix):
                list_of_files[filename] = os.sep.join([dirpath, filename])


output_path = os.sep.join([path,'sum_manta_result.txt'])
output = open(output_path,'w')

header_list=["SAMPLE","CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","FORMAT_DATA"]
output.write('\t'.join(header_list)+'\n')

for file in list_of_files.values():
    with gzip.open(file,'rt') as fp:

        for line in fp:
            if line[0] == '#' and line[1] != '#':
                samplename = line.split()[9]
            elif line[0] != '#':
                new_line = samplename + "\t" + line
                output.writelines(new_line)

output.close()