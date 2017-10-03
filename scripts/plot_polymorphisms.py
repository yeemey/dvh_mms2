#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 13:15:18 2017

@author: ymseah
"""

from lost_polymorphisms import ComparePolymorphisms
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

input_directory = '/opt/data/wkim-rsrch/breseq_results/'
output_directory = '/home/NETID/ymseah/Projects/plot_polymorphisms/results/'
ancestor_gd_path = input_directory + 'sic_Ancestor_breseq/output/0.gd'
evolution_lines = ['HA3', 'HE3', 'HR2', 'HS3', 'UA3', 'UE3', 'UR1', 'US1']

def subset_monotonic(df):
    '''
    Returns 3 lists of column names that have monotonically increasing, decreasing, and non-monotonic values. 
    '''
    mono_inc = []
    mono_dec = []
    mono_oth = []
    for col in df.columns:
        if df[col].is_monotonic_increasing:
            mono_inc.append(col)
        elif df[col].is_monotonic_decreasing:
            mono_dec.append(col)
        else:
            mono_oth.append(col)
    return mono_inc, mono_dec, mono_oth

def plot_line_summary(line, df, plottitle):
    data = df.loc[line.upper()].dropna(axis=1, how='all').fillna(0)
    inc, dec, oth = subset_monotonic(data)
    fig, axes = plt.subplots(3,1, sharex=True)
    if len(inc) > 0:
        data[inc].plot(ax=axes[0], legend=False, title=plottitle)
    else:
        axes[0].set_title(plottitle)
    if len(dec) > 0:
        data[dec].plot(ax=axes[1], legend=False)
    if len(oth) > 0:
        data[oth].plot(ax=axes[2], legend=False)
    plt.savefig(plottitle + '_summary.png')

def plot_line_individual(line, df, genome_id):
    data = df.loc[line.upper()].dropna(axis=1, how='all').fillna(0)
    inc, dec, oth = subset_monotonic(data)
    if len(inc) > 0:
        data[inc].plot(title = line.upper() + ':' + genome_id + 'is_monotonic_increasing', figsize=(20,15)).legend(loc='upper left', ncol=6)
        plt.savefig(line.upper() + '_' + genome_id + '_monotonic_increasing.png')
    if len(dec) > 0:
        data[dec].plot(title = line.upper() + ':' + genome_id + 'is_monotonic_decreasing', figsize=(20,15)).legend(loc='upper right', ncol=6)
        plt.savefig(line.upper() + '_' + genome_id + '_monotonic_decreasing.png')
    if len(oth) > 0:
        data[oth].plot(title = line.upper() + ':' + genome_id + ' non-monotonic', figsize=(20,15)).legend(loc='upper right', ncol=6)
        plt.savefig(line.upper() + '_' + genome_id + '_non-monotonic.png')

'''
# frequencies from HTML output
compare_files = glob.glob(input_directory + 'compare/[H|U]*.html')
print('All gdtools COMPARE html files: ')
print(compare_files)

cp = ComparePolymorphisms()
for file in compare_files:
    print(file)
    evolution_line = file[-8:-5]
    print(evolution_line)
    bs_from_html = cp.parse_compare_html(file)
    mutation_freqs_dict = cp.get_html_generation_frequencies(bs_from_html)
    output_csv = output_directory + evolution_line + '_html_freqs.csv'
    html_freqs_df = cp.html_frequencies_to_df(mutation_freqs_dict, save_csv=True, csv_filename=output_csv)
'''

# frequencies from GD output
cp_gd = ComparePolymorphisms()
for line in evolution_lines:
    line_df = cp_gd.get_all_gd(line, input_directory, ancestor_gd_path)
    output_csv = output_directory + line + '_gd_freqs.csv'
    frequencies_df = cp_gd.gd_frequencies_to_df(line_df, save_csv=True, csv_filename=output_csv)

# gather all frequencies (csv) into multi-indexed dataframe
all_csv = glob.glob('/Users/ymseah/Repositories/dataviz/data/*_gd_freqs.csv')
df_list = []
for csv in all_csv:
    df = pd.read_csv(csv)
    df_list.append(df)
all_df = pd.concat(df_list, ignore_index=True)
all_df['polymorphism'] = all_df['position'].astype(str).str.cat(all_df['entry_type'], sep=' ').str.cat(all_df['mutation_detail'], sep=' ')
all_df = all_df[['line', 'genome_id', 'generation', 'polymorphism', 'polymorphism_frequency']].copy()
all_df_pvt = all_df.pivot_table(index=['genome_id', 'line', 'generation'], 
                                columns='polymorphism', values='polymorphism_frequency')
# dataframes by organism
dvh_df = all_df_pvt.loc['NC_002937'].dropna(axis=1, how='all')
mm_df = all_df_pvt.loc['NC_005791'].dropna(axis=1, how='all')
dvplasmid_df = all_df_pvt.loc['NC_005863'].dropna(axis=1, how='all')

# plot polymorphisms by monotonic/non-monotonic frequency trajectories
for line in all_df_pvt.index.levels[1]:
    plot_line_summary(line, dvh_df, (line + '-NC_002937'))
    plt.close()
    plot_line_individual(line, dvh_df, 'NC_002937')
    plt.close()
    plot_line_summary(line, mm_df, (line + '-NC_005791'))
    plt.close()
    plot_line_individual(line, mm_df, 'NC_005791')
    plt.close()
    plot_line_summary(line, dvplasmid_df, (line + '-NC_005863'))
    plt.close()
    plot_line_individual(line, dvplasmid_df, 'NC_005863')
    plt.close()
    