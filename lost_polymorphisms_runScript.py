#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 14:28:40 2017

@author: ymseah
"""
from lost_polymorphisms import ComparePolymorphisms
import glob

input_directory = '/opt/data/wkim-rsrch/breseq_results/'
output_directory = '/home/NETID/ymseah/Projects/lost_polymorphisms_in_breseq/results/'
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
    mutation_freqs_dict = cp.get_generation_frequencies(bs_from_html)
    suspect_freqs_dict = cp.get_suspect_frequencies(mutation_freqs_dict)
    annotated_gd_files = glob.glob(input_directory + 'sic_' + evolution_line + '*/output/*.gd')
    print(annotated_gd_files)
    all_dataframes = [ancestor_df]
    for genome_diff in annotated_gd_files:
        print(genome_diff)
        generation = int(genome_diff[-8:-3].split('-')[1])
        print(generation)
        dataframe = cp.annotated_gd_to_df(genome_diff, generation)
        all_dataframes.append(dataframe)
    print(all_dataframes)
    evolution_line_dataframe = cp.summary_df(evolution_line, all_dataframes, output_directory)
    new_suspect_freqs_dict = cp.get_reject_reasons(evolution_line_dataframe, suspect_freqs_dict)
    cp.write_frequency_dicts_to_file(suspect_freqs_dict, output_directory + evolution_line + '_suspect')
    cp.write_frequency_dicts_to_file(new_suspect_freqs_dict, output_directory + evolution_line + '_suspect_reasons')     