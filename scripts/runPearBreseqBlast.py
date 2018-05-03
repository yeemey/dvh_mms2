#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 16:48:53 2016

python 3.5.2
BioPython 1.68

@author: ymseah
"""

from blastUnmapped_svr import run_software, parse_results
import re, os, subprocess

'''
sample names and directory paths
'''
sample_name = 'sic_UA3_S2_L001'
forward_read = sample_name + '_R1_001.fastq'
reverse_read = sample_name + '_R2_001.fastq'
scy_sic_output = '/home/NETID/ymseah/Projects/20180307_adapters_2015/results/scy_sic/'
pear_arg = '~/./pear'
pear_output = '/home/NETID/ymseah/Projects/20180307_adapters_2015/results/pear/'
breseq_arg = '~/breseq'
ref_dir = ''
breseq_output = '/home/NETID/ymseah/Projects/20180307_adapters_2015/results/breseq/'
gdtools_arg = '~/gdtools'
batch_reads_dir = '/opt/data/wkim-data/FilteredFastQC/'

'''
run 
'''
run_object = run_software()
#print(run_object.get_all_pear_samples('/Users/ymseah/scriptingProjects/python/test/'))
#run_object.batch_run_pear(pear_arg, batch_reads_dir, pear_output)
#run_object.batch_run_breseq(breseq_arg, ref_dir, pear_output, breseq_output)
run_object.run_pear(pear_arg, scy_sic_output + forward_read, scy_sic_output + reverse_read, pear_output, sample_name)
run_object.run_breseq(breseq_arg, ref_dir, pear_output, sample_name, breseq_output)
