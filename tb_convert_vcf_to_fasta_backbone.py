#!/usr/bin/env python

# standard library imports
import os
import sys
import argparse
from pathlib import Path
import json
import gzip

parser = argparse.ArgumentParser(description='ADD YOUR DESCRIPTION HERE')
parser.add_argument('-i', '--input', help='Input folder name', required=True)
parser.add_argument('-o', '--output', help='Absolute path of Output file', required=True)
parser.add_argument('-t', '--target', help='Extension of target input file Ex. ".vcf", ".vcf.gz"', required=True)
parser.add_argument('-s', '--start_line', help='the first three characters for indication where a data line starts.',
                    required=False)

args = parser.parse_args()


def check_total_samples(sample_list):
    """Basie function.
    This function is for checking the total number of sample in the input list.

    Args:
        sample_list (list): The list of samples.

    Return:
        int: The total sample number in the list.

    """
    return len(sample_list)


def get_in_file_list(in_dir, file_name_part, path_num=1, alternative_directory=None):
    """Basie function.
    This function is for getting a list of files in the target location based on part of a file name.

    Args:
        in_dir (string): a name of a target folder.
        file_name_part (string): a end-part of the file name that you want to get a list on.
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


def prep_output(out_dir, out_file_name=None, path_num=1, alternative_directory=None):
    """Basic function.
    This function is for preparing an output at a target location folder and
        also creating the output folder if the target does not already exist.

    Args:
        out_dir (string): A name of a output folder
        out_file_name (string): A name of a output file
        path_num (int): option for using alternative output folder part. 1 is default.
        alternative_directory (string): if you want to output a result to a different location
                                        not the same location with this script, you can input
                                        the path with this parameter.

    Return:
        string: a path to output folder.

    """
    if path_num == 1:
        cwd = os.path.dirname(os.path.realpath(sys.argv[0]))
    else:
        cwd = alternative_directory
    out_path = os.path.join(cwd, out_dir)
    if not os.path.exists(out_path):
        os.system("mkdir " + out_path)
    if out_file_name == None:
        return out_path;
    else:
        out_file_path = os.path.join(out_path, out_file_name)
        return out_file_path;


def extra_name(some_string, codition=1):
    """Basie function.
    This function extracts a file name from a file path in different conditions.

    Args:
        some_string (string): A input file path.
        condition (int): The condition for extracting the file name from the path.
                        1 - a string between '/' and '.' (default)
                        2 - a string between '/' and '_'
                        3 - a string between start to '.'
                        4 - a string between start to '_'

    Return:
        string: a string name.

    """
    first_position = 0
    last_position = 0
    unders_position = 0
    size = len(some_string)
    for i in range(size):
        if some_string[(size - 1) - i] == '/' and first_position == 0:
            first_position = (size - i)
        if some_string[(size - 1) - i] == '_' and unders_position == 0:
            unders_position = (size - i - 1)
        if some_string[i] == '.' and last_position == 0:
            last_position = i
    if codition == 1:
        return some_string[first_position:last_position];
    elif codition == 2:
        return some_string[first_position:unders_position];
    elif codition == 3:
        return some_string[:last_position];
    elif codition == 4:
        return some_string[:unders_position];


def create_key_position_list(in_dir, input_file_extension, sample_start_line='NC_', out_file='Master_list',
                             out_file_name1='master_aon_list.json', out_file_name2='number_SNPs_per_samples.txt', output_file_path=''):
    """Advanced function.
    This function is for getting a list of all SNP position from all the isolates
    and dropping exclude positions out.

    To use this function:
        - Check in the input file for first three characters for indicating a result data line.
        - This function is only worked on result file from SAMtools program (.vcf).
        - The exclude position file much be a position number per a line.
        - This function will drop SNPs that have more than one types of mutation.
        - This function only output all position of SNPs in a population for future reference.
        - All input files must be in the input folder.
        - Another output is a list that shows total SNPs from each isolate.

    Args:
        in_dir (string): a name of a input target folder.
        input_file_extension (string): file extension of input file.
        in_ex_posi (string): a name of a exclude position list file.
        out_file (string): a name of SNP output position list file.
        out2_file (string): a total number SNP for each isolate file name.
        sample_start_line (string): the first three characters for indication where a data line starts.

    Return:
        None
    """
    f_list = get_in_file_list(in_dir, input_file_extension)
    outpath1 = prep_output(out_file, out_file_name1, 2, output_file_path)
    print(outpath1)
    outpath2 = prep_output(out_file, out_file_name2, 2, output_file_path)

    ms_list = {}
    # count the total of the SNPs from all the samples
    mas_cout = 0

    with open(outpath2, 'w') as o2:
        for j in f_list:
            all_snps = 0
            print(j)
            if j.endswith(".gz"):
                with gzip.open(j, 'rt') as f:
                    for line in f:
                        if line[0:3] == sample_start_line:
                            line = line.split(';')
                            line2 = line[0].split("\t")
                            all_snps += 1
                            ms_key = ms_list.keys()
                            if (int(line2[1]) not in ms_key):
                                mas_cout += 1
                                ms_list[int(line2[1])] = str(line2[3])
                    p = 'GZ Number of SNPs in ' + extra_name(j, 1) + ' : ' + str(all_snps) + '\n'
                    o2.write(p)
            elif j.endswith(".vcf"):
                with open(j, 'r') as f:
                    for line in f:
                        if line[0:3] == sample_start_line:
                            line = line.split(';')
                            line2 = line[0].split("\t")
                            all_snps += 1
                            ms_key = ms_list.keys()
                            if (int(line2[1]) not in ms_key):
                                mas_cout += 1
                                ms_list[int(line2[1])] = str(line2[3])
                    p = 'Number of SNPs in ' + extra_name(j, 1) + ' : ' + str(all_snps) + '\n'
                    o2.write(p)
        o2.write('Total number of SNP list ' + str(len(ms_list)))
    ms_key = ms_list.keys()
    ms_key = sorted(ms_key)
    with open(outpath1, 'w') as y:
        json.dump(ms_list, y, sort_keys=True, indent=1)     ## save master template for later use
        #for h in ms_key:
            #y.write(str(h) + '\n')
    return ms_list;


def create_SNPs_seq(master_SNPs, in_dir, input_file_extension,
                    output_file_name, sample_start_line, output_file_path):
    """Advanced function.
    This function takes in several input files to build SNPs sequence in a Fasta format.

    To use this function:
        - Check in the input file for first three characters for indicating a result data line.
        - This function is only worked on result file from SAMtools program (.vcf).
        - The exclude position file much be a position number per a line.
        - This function will drop SNPs that have more than one types of mutation.
        - All input files must be in the input folder.
        - The whole genome sequence file name must be the same with VCF file name.
        - The output file will be in fasta file format.
        - The function will add pre-assigned group extension at the end of isolate name.

    Args:
        in_vcf_dir (string): a name of a SNPs data input folder.
        in_vcf_file_extension (string): file extension of SNPs data file.
        in_WGenome_dir (string): a name of a artificial whole genome sequence folder.
        in_WGenome_extension (string): file extension of whole genome sequence file.
        master_SNPs_position_file (string): a name of all SNPs position list file.
        isolate_name_file (string): a name of isolate pre-assign group file
        in_ex_posi (string): a name of a exclude position list file.
        output_file_name (string): a name of SNPs seq output file.
        sample_start_line (string): the first three characters for indication where a data line starts.

    Return:
        None
    """
    first_time_flag = True
    #Master snp template output path
    master_snp_file = prep_output(out_dir='', out_file_name='master_snp_template.json', path_num=2, alternative_directory=output_file_path)

    # output path
    out_file = prep_output(out_dir='', out_file_name=output_file_name, path_num=2, alternative_directory=output_file_path)
    # master list file
    ms_key = master_SNPs.keys()
    ms_key = sorted(ms_key)
    # get input file list
    f_list = get_in_file_list(in_dir, input_file_extension)
    # count the number of sample for checking the progress of the program
    count_total = 0
    # get total number of samples
    total = check_total_samples(f_list)
    # process to create SNPs seq from all the components
    for file_snp in f_list:
        # show the progace fo the program and the name of the files the program work on
        count_total += 1
        print('------------sample mumber : ' + str(count_total) + ' from ' + str(total) + '-------------------')
        print(file_snp)
        # create the temp seq dictionary by position number of SNPs and fill with '-'
        sample_seq_dic = {}
        for p in ms_key:
            sample_seq_dic[int(p)] = '-'
        print(len(sample_seq_dic))
        print('Dic_done')
        # create SNPs list [SNP_positon , SNP_base]
        ms_list = []
        if file_snp.endswith(".gz"):
            with gzip.open(file_snp, 'rt') as f:
                count = 0
                for line in f:
                    if line[0:3] == sample_start_line:
                        count += 1
                        line = line.split(';')
                        line2 = line[0].split("\t")
                        # print line2
                        ms_list.append([int(line2[1]), str(line2[4])])
        elif file_snp.endswith(".vcf"):
            with open(file_snp, 'r') as f:
                count = 0
                for line in f:
                    if line[0:3] == sample_start_line:
                        count += 1
                        line = line.split(';')
                        line2 = line[0].split("\t")
                        # print line2
                        ms_list.append([int(line2[1]), str(line2[4])])
        print(len(ms_list))
        print('Snp_done')
        # fill temp seq dictionary with SNPs by position
        for b in ms_list:
            sample_seq_dic[b[0]] = b[1]
        print('SNPs_fill_done')
        # fill temp seq dictionary with base by position
        for i in ms_key:
            if (sample_seq_dic[i] == '-'):
                sample_seq_dic[i] = master_SNPs[i]
        print('Ref_SNPs_fill_done')
        # create empty list for put SNPS seq together
        seq_snp_sample = []
        # check nember of seq bafore and after put them together
        for l in ms_key:
            seq_snp_sample.append(str(sample_seq_dic[l]))
        print('len seq before join')
        print(len(seq_snp_sample))
        seq_snp_sample = ''.join(seq_snp_sample)
        print('len seq after join')
        print(len(seq_snp_sample))
        # Put SNPs seq them in the output file
        if first_time_flag == True:
            with open(out_file, 'w') as f:
                f.write('>' + str(extra_name(file_snp, 1)) + '\n')
                f.write(seq_snp_sample + '\n')
                first_time_flag = False
        else:
            with open(out_file, 'a') as f:
                f.write('>' + str(extra_name(file_snp, 1)) + '\n')
                f.write(seq_snp_sample + '\n')
    return;


# --------main---------
if __name__ == "__main__":
    if args.start_line == None:
        #start_line = 'NC_'
        start_line = 'Chr'
    else:
        start_line = args.start_line

    output_absolute_path = Path(args.output)
    output_filename = output_absolute_path.name
    output_dir = output_absolute_path.parent

    output_dir.mkdir(parents=True, exist_ok=True)    ## create output dir that user specified if not exist

    mast_position_list = create_key_position_list(in_dir=args.input, input_file_extension=args.target,
                                                  sample_start_line=start_line,
                                                  out_file='Master_list', out_file_name1='master_list.json',
                                                  out_file_name2='number_SNPs_per_samples.txt', output_file_path=output_dir)
    create_SNPs_seq(mast_position_list, in_dir=args.input,
                    input_file_extension=args.target,
                    output_file_name=output_filename,
                    sample_start_line=start_line,
                    output_file_path=output_dir)
