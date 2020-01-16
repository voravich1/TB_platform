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
path="/Volumes/10TBSeagateBackupPlus/NBT/thalassemia/analysis/HS08/analysis/MantaSVvcf"


list_of_files = {}
for (dirpath, dirnames, filenames) in os.walk(path):
    for filename in filenames:
        if filename.endswith('_filtered.vcf.gz'):
            list_of_files[filename] = os.sep.join([dirpath, filename])


output_path = os.sep.join([path,'sum_manta_result.txt'])
output = open(output_path,'w')


for file in list_of_files.values():
    with gzip.open(file,'rt') as fp:

        for line in fp:
            if line[0] == '#' and line[1] != '#':
                samplename = line.split()[9]
            elif line[0] != '#':
                new_line = samplename + "\t" + line
                output.writelines(new_line)

output.close()