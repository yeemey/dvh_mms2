#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 11:22:47 2017

@author: ymseah
"""

from lost_polymorphisms import ComparePolymorphisms
import glob

evolution_line = input('Evolution line (e.g., HA3): ')
input_directory = '/opt/data/wkim-rsrch/breseq_results/'
compare_input = input_directory + 'compare/' + evolution_line + '.html'
output_directory = '/home/NETID/ymseah/Projects/lost_polymorphisms_in_breseq/results/'
generations = [0, 100, 300, 500, 780, 1000]

print('COMPARE html file: ' + compare_input)

cp = ComparePolymorphisms()
ancestor_df = cp.annotated_gd_to_df(input_directory + 'sic_Ancestor_breseq/output/0.gd', 0)

bs_from_html = cp.parse_compare_html(compare_input)
mutation_freqs_dict = cp.get_html_generation_frequencies(bs_from_html)
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
cp.write_html_frequency_dicts_to_file(mutation_freqs_dict, output_directory + evolution_line)
evidence_dict = cp.get_polymorphism_evidence(evolution_line_dataframe, mutation_freqs_dict)
cp.write_evidence_dicts_to_file(evidence_dict, output_directory + evolution_line)