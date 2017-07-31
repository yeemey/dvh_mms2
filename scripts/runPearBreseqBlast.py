#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:48:53 2016

python 3.5.2
BioPython 1.68

@author: ymseah
"""

from blastUnmapped import run_software, parse_results
import re, os, subprocess

run_object = run_software()
input_dir = '/home/NETID/ymseah/Projects/ngs_preprocess/20170619/results/'
forward_read = 'sic_UA3_S2_L001_R1_001.fastq'
reverse_read = 'sic_UA3_S2_L001_R2_001.fastq'
sample_name = 'sic_UA3_S2_L001'
output_dir = '/home/NETID/ymseah/Projects/ngs_preprocess/20170619/results/'
batch_reads_dir = '/opt/data/wkim-data/FilteredFastQC/'
#print(run_object.get_all_pear_samples('/Users/ymseah/scriptingProjects/python/test/'))
#run_object.batch_run_pear('/opt/data/wkim-data/FilteredFastQC/'', 
#                          'home/NETID/ymseah/Projects/Low_Mapping_in_breseq/results/20161201/')
#run_object.batch_run_breseq('/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/results/20161202/', 
#                            '/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/results/breseq_results/')
#run_object.run_breseq('/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/results/20161202/',
#                      'sic_HE3.118','/home/NETID/ymseah/Projects/Low_Mapping_in_breseq/results/breseq_results/')
#run_object.run_gdtools_compare('sic_HE3-15', 'sic_HE3.45', 'sic_HE3-76', 'sic_HE3-118')
#run_object.run_blastn_remote('sic_HA3.45')
run_object.run_pear(input_dir + forward_read, input_dir + reverse_read, output_dir, sample_name)
run_object.run_breseq(input_dir, sample_name, output_dir)

