import json
import vcfpy

sv_db_file = "/Users/worawich/Download_dataset/TB_platform_test/test_data/test_data_bgi/lin_db/sv_lin_db.txt"

with open(sv_db_file,'r') as sv_db_reader:
    for line in sv_db_reader:
