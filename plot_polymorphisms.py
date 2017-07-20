#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 13:15:18 2017

@author: ymseah
"""

from lost_polymorphisms import ComparePolymorphisms
import glob

input_directory = '/opt/data/wkim-rsrch/breseq_results/'
output_directory = '/home/NETID/ymseah/Projects/plot_polymorphisms/results/'
generations = [0, 100, 300, 500, 780, 1000]

compare_files = glob.glob(input_directory + 'compare/[H|U]*.html')
print('All gdtools COMPARE html files: ')
print(compare_files)

cp = ComparePolymorphisms()
ancestor_df = cp.annotated_gd_to_df(input_directory + 'sic_Ancestor_breseq/output/0.gd', 0)
for file in compare_files:
    print(file)
    evolution_line = file[-8:-5]
    print(evolution_line)
    bs_from_html = cp.parse_compare_html(file)
    mutation_freqs_dict = cp.get_html_generation_frequencies(bs_from_html)
    output_csv = output_directory + evolution_line + '_html_freqs.csv'
    html_freqs_df = cp.html_frequencies_to_df(mutation_freqs_dict, save_csv=True, csv_filename=output_csv)
    
