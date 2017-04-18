#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 15:39:55 2017

@author: ymseah
"""

from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
import re

#Frequency data across generations
ua3 = BeautifulSoup(open('compare/UA3.html'), 'html.parser')
ue3 = BeautifulSoup(open('compare/UE3.html'), 'html.parser')
ur1 = BeautifulSoup(open('compare/UR1.html'), 'html.parser')
us1 = BeautifulSoup(open('compare/US1.html'), 'html.parser')
ha3 = BeautifulSoup(open('compare/HA3.html'), 'html.parser')
he3 = BeautifulSoup(open('compare/HE3.html'), 'html.parser')
hr2 = BeautifulSoup(open('compare/HR2.html'), 'html.parser')
hs3 = BeautifulSoup(open('compare/HS3.html'), 'html.parser')

#1. Identify ref genome and polymorphism position
#2. Identify frequencies across generations
#3. Link 1. and 2.
def parse_compare_html(filepath):
    pass

def get_generation_frequencies(data_html_object):
    '''
    Takes compare html object parsed by BeautifulSoup, outputs dictionary of mutation frequencies across generations
    Key: (reference genome ID, position of mutation)
    Value: [mutation, ancestor frequency, generation 100 f, gen 300 f, gen 500 f, gen 780 f, gen 1000 f]
    '''
    table_cells = data_html_object.body.find_all('td')
    generation_frequencies_dict = {}
    for cell in table_cells:
        if re.match('NC_', str(cell.string)):
            position = cell.next_sibling.next_sibling.next_sibling
            clean_position = int(re.sub(',', '', position.string))
            gen_freqs_key = (cell.string, clean_position)
            mutation = position.next_sibling.next_sibling.next_sibling
            freq_ancestor = mutation.next_sibling.next_sibling.next_sibling
            freq_gen100 = freq_ancestor.next_sibling
            freq_gen300 = freq_gen100.next_sibling
            freq_gen500 = freq_gen300.next_sibling
            freq_gen780 = freq_gen500.next_sibling
            freq_gen1000 = freq_gen780.next_sibling
            gen_freqs_value = [mutation.string, freq_ancestor.string, freq_gen100.string, freq_gen300.string, freq_gen500.string, freq_gen780.string, freq_gen1000.string]
            generation_frequencies_dict[gen_freqs_key] = gen_freqs_value
    return generation_frequencies_dict

#5. Identify a-n-a pattern.

suspect_frequencies_dict = {}
with open('suspect_frequencies.tsv', 'w') as output_file:
    output_file.write('ref_genome\tposition\tmutation\tfreq_anc\tfreq_100\tfreq_300\tfreq_500\tfreq_780\tfreq_1000\n')
    for key, value in generation_frequencies_dict.items():
        counter = 2
        while counter <= len(value):
            if counter + 1 < len(value):
                if value[counter] == None:
                    if value[counter - 1] == '100%' or value[counter + 1] == '100%':
                        suspect_frequencies_dict[key] = value
                        output_file.write(key[0] + '\t' + str(key[1]) + '\t')
                        for item in value:
                            output_file.write(str(item) + '\t')
                        output_file.write('\n')
                elif value[counter] == '100%':
                    if value[counter - 1] == None or value[counter + 1] == None:
                        suspect_frequencies_dict[key] = value
                        output_file.write(key[0] + '\t' + str(key[1]) + '\t')
                        for item in value:
                            output_file.write(str(item) + '\t')
                        output_file.write('\n')
            counter += 2

UA3_100 = pd.read_table('/Users/ymseah/Documents/sic_UA3-15_breseq/annotated.gd', comment='#', names=range(50), dtype=str)
UA3_100.insert(0, 'generation', 100)
UA3_300 = pd.read_table('/Users/ymseah/Documents/sic_UA3.45_breseq/annotated.gd', comment='#', names=range(50), dtype=str)
UA3_300.insert(0, 'generation', 300)
UA3_500 = pd.read_table('/Users/ymseah/Documents/sic_UA3-76_breseq/annotated.gd', comment='#', names=range(50), dtype=str)
UA3_500.insert(0, 'generation', 500)
UA3_780 = pd.read_table('/Users/ymseah/Documents/sic_UA3.118_breseq/annotated.gd', comment='#', names=range(50), dtype=str)
UA3_780.insert(0, 'generation', 780)
UA3_1000 = pd.read_table('/Users/ymseah/Documents/sic_UA3_S2_L001_breseq/output/evidence/annotated.gd', comment='#', names=range(50), dtype=str)
UA3_1000.insert(0, 'generation', 1000)

UA3_df  = pd.concat([UA3_100, UA3_300, UA3_500, UA3_780, UA3_1000], ignore_index=True)
UA3_df.insert(0, 'line', 'UA3')
UA3_df.insert(2, 'frequency', 0.0)
UA3_df.insert(3, 'gene_product', '')
UA3_df.insert(4, 'gene_position', '')
UA3_df.insert(5, 'reject', '')

for row in UA3_df.itertuples():
    #check each column
    col_index = 6
    while col_index < 50:
        #1. polymorphism frequencies
        if re.match('frequency=', str(UA3_df.loc[row[0], col_index])):
            UA3_df.loc[row[0], 'frequency'] = re.sub('frequency=', '', str(UA3_df.loc[row[0], col_index]))
        #2. gene products
        elif re.match('gene_product=', str(UA3_df.loc[row[0], col_index])):
            UA3_df.loc[row[0], 'gene_product'] = re.sub('gene_product=', '', str(UA3_df.loc[row[0], col_index]))
        #3. polymorphism rejection reasons
        elif re.match('reject=', str(UA3_df.loc[row[0], col_index])):
            UA3_df.loc[row[0], 'reject'] = re.sub('reject=', '', str(UA3_df.loc[row[0], col_index]))
        #4. gene annotations
        elif re.match('gene_position=', str(UA3_df.loc[row[0], col_index])):
            UA3_df.loc[row[0], 'gene_position'] = re.sub('gene_position=', '', str(UA3_df.loc[row[0], col_index]))
        col_index += 1
    #set frequencies type to float
    if re.match('1|2|3|4|5|6|7|8|9', str(UA3_df.loc[row[0], 'frequency'])):
        UA3_df.loc[row[0], 'frequency'] = float(UA3_df.loc[row[0], 'frequency'])
    else:
        UA3_df.loc[row[0], 'frequency'] = 0.0
    #set positions (col 4) type to int
    UA3_df.loc[row[0], 4] = int(UA3_df.loc[row[0], 4])
    #set reject col to 'NA' when no reject reason given.
    if (UA3_df.loc[row[0], 'reject'] == '') & (UA3_df.loc[row[0], 2] == '.'):
        UA3_df.loc[row[0], 'reject'] = 'NA'

UA3_df.rename(columns = {0: 'entry_type', 1: 'item_id', 2: 'evidence_ids', 3: 'ref_genome', 4:'position'}, inplace=True)
ua3df_subset = UA3_df[['line', 'generation', 'frequency', 'gene_product', 'gene_position', 'reject', 'entry_type', 'item_id', 'evidence_ids', 'ref_genome', 'position']].copy()
ua3df_subset.to_csv('/Users/ymseah/Documents/ua3.csv', index=False)

#6: Identify reasons for a-n-a frequencies.

ua3df_subset_mutations_only = ua3df_subset[(ua3df_subset['entry_type'] == 'SNP') | 
                                           (ua3df_subset['entry_type'] == 'SUB') | 
                                           (ua3df_subset['entry_type'] == 'DEL') | 
                                           (ua3df_subset['entry_type'] == 'INS') | 
                                           (ua3df_subset['entry_type'] == 'MOB') | 
                                           (ua3df_subset['entry_type'] == 'AMP') | 
                                           (ua3df_subset['entry_type'] == 'CON') | 
                                           (ua3df_subset['entry_type'] == 'INV')]
for key, value in suspect_frequencies_dict.items():
    row_indices = ua3df_subset_mutations_only[(ua3df_subset_mutations_only['ref_genome'] == key[0]) & 
                                              (ua3df_subset_mutations_only['position'] == key[1])].index.tolist()
    print(key, row_indices)
    for row in row_indices:
        reject_reason = ua3df_subset_mutations_only.loc[row, 'reject']
        generation = ua3df_subset_mutations_only.loc[row, 'generation']
        entry_type = ua3df_subset_mutations_only.loc[row, 'entry_type']
        print(row, reject_reason, generation, entry_type)
        if reject_reason != '':
            if generation == 100:
                value[2] = str(value[2]) + ' ' + reject_reason
            elif generation == 300:
                value[3] = str(value[3]) + ' ' + reject_reason
            elif generation == 500:
                value[4] = str(value[4]) + ' ' + reject_reason
            elif generation == 780:
                value[5] = str(value[5]) + ' ' + reject_reason
            elif generation == 1000:
                value[6] = str(value[6]) + ' ' + reject_reason

for key, value in suspect_frequencies_dict.items():
    print(key, value)
