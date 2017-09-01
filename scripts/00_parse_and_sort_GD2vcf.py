#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 13:04:05 2017

@author: ymseah
"""
#from sys import argv
#gd_file, sample_name = argv
import re
gd_file = '/Users/ymseah/Documents/sic_UA3-15_breseq/annotated.gd'

def create_mock_vcf(vcf_in):
    with open('mock.vcf', 'w') as outfile:
        outfile.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t\n')
        with open(vcf_in, 'r') as infile_vcf:
            for line in infile_vcf:
                if line.startswith('##') or line.startswith('#'):
                    continue                
                outfile.write(line.strip('\n') + '\tAD:DP\t\n')
    return            

def get_gd_read_counts(string):
    line = string.split('\t')
    for data in line:
        if re.match('reject', data):
            return '.'
        else:
            continue
    line_data = [data.strip() for data in line if re.search('=', data)]
    line_key_data_split = [item.split('=') for item in line_data]
    line_key_data_dict = {}
    for each in line_key_data_split:
        line_key_data_dict[each[0]] = each[1]
    
    if line[0] == 'RA':
        major_cov_per_strand = line_key_data_dict['major_cov']
        major_cov = major_cov_per_strand.split('/')
        total_major_cov = int(major_cov[0]) + int(major_cov[1])
        minor_cov_per_strand = line_key_data_dict['minor_cov']
        minor_cov = minor_cov_per_strand.split('/')
        total_minor_cov = int(minor_cov[0]) + int(minor_cov[1])
        AD = str(total_major_cov) + ',' + str(total_minor_cov)
        DP = str(total_major_cov + total_minor_cov)    
    elif line[0] == 'JC':
        AD = line_key_data_dict['new_junction_read_count']
        DP = AD
    elif line[0] == 'MC':
        left_outside_cov = line_key_data_dict['left_outside_cov']
        right_outside_cov = line_key_data_dict['right_outside_cov']
        if left_outside_cov != 'NA' and right_outside_cov != 'NA':
            AD = int(left_outside_cov) + int(right_outside_cov)
        elif left_outside_cov != 'NA' and right_outside_cov == 'NA':
            AD = int(left_outside_cov)
        elif left_outside_cov == 'NA' and right_outside_cov != 'NA':
            AD = int(right_outside_cov)
        else:
            AD = 0
        DP = AD
    
    AD_DP_counts = str(AD) + ':' + str(DP)
    return AD_DP_counts

def add_counts_to_vcf(string, counts_dict):
    line = string.split('\t')
    if (line[0], line[1]) in counts_dict:
        string.strip('\n') + counts_dict[line[0], line[1]] + '\n'
    return string
    
create_mock_vcf('/Users/ymseah/Repositories/rice_desai/sic_UA3-15.vcf')
gd_read_depths = {}
with open(gd_file, 'r') as infile_gd:    
    for gd_line in infile_gd:
        if gd_line.startswith('RA') or gd_line.startswith('MC'):
            read_depth = get_gd_read_counts(gd_line)
            gd_line_data = gd_line.split('\t')
            gd_read_depths[(gd_line_data[3], gd_line_data[4])] = read_depth
        elif gd_line.startswith('JC'):
            
print(gd_read_depths)