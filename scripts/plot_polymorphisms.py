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

def fluctuating_bools(df):
    '''
    Returns 2 lists of alternating bool values. Default list length is the number of indexes in df.
    '''
    true_false = []
    false_true = []
    count = 0
    while count < len(df.index.tolist()):
        true_false.append(True)
        false_true.append(False)
        count += 1
        if count < len(df.index.tolist()):
            true_false.append(False)
            false_true.append(True)
            count += 1
    return true_false, false_true

def greater_than_last_gen(df):
    '''Compares if frequencies in each row greater than in previous row, 
    returns df of bool values'''
    comparisons = {}
    count = 0
    while count < len(df.index)-1:
        comparisons[str(count+1) + '_to_' + str(count)] = df.loc[df.index[count+1]] > df.loc[df.index[count]]
        count += 1
    compare_df = pd.DataFrame(comparisons).T
    return compare_df

def subset_fluctuating_frequencies(df):
    '''Returns list of df columns with fluctuating polymorphisms'''
    t_f, f_t = fluctuating_bools(df)
    fluctuating_polymorphisms = []
    for col in df.columns:
        if (df[col].tolist() == t_f) | (df[col].tolist() == f_t):
            fluctuating_polymorphisms.append(col)
    return fluctuating_polymorphisms

def subset_some_fluctuation(df, tf_len=3):
    t_f, f_t = fluctuating_bools(df)
    some_t_f = t_f[:tf_len]
    some_f_t = f_t[:tf_len]
    some_flux = []
    for col in df.columns:
        bools_list = df[col].tolist()
        count = 0
        while count <= (len(bools_list)-3):
            bools_slice = bools_list[count:(count+3)]
            if (bools_slice == some_t_f) | (bools_slice == some_f_t):
                some_flux.append(col)
            count += 1
    return list(set(some_flux))

def subset_monotonic(df):
    '''
    Returns 2 lists of column names that have monotonically increasing, and decreasing values. 
    '''
    mono_inc = []
    mono_dec = []
    for col in df.columns:
        if df[col].is_monotonic_increasing:
            mono_inc.append(col)
        elif df[col].is_monotonic_decreasing:
            mono_dec.append(col)
    return mono_inc, mono_dec

# Plotting functions

def plot_fluctuating(df, some_flux=False, flux_len=3):
    compare_freqs = greater_than_last_gen(df)
    if some_flux == False:
        flux = subset_fluctuating_frequencies(compare_freqs)
    else:
        flux = subset_some_fluctuation(compare_freqs, tf_len = flux_len)
    df[flux].plot().legend(loc='upper right')   

def plot_line_summary(line, df, plottitle):
    data = df.loc[line.upper()].dropna(axis=1, how='all').fillna(0)
    # monotonic
    inc, dec = subset_monotonic(data)
    fig, axes = plt.subplots(5,1, sharex=True)
    if len(inc) > 0:
        data[inc].plot(ax=axes[0], legend=False, title=plottitle)
    else:
        axes[0].set_title(plottitle)
    if len(dec) > 0:
        data[dec].plot(ax=axes[1], legend=False)
    # all fluctuating
    compare_freqs = greater_than_last_gen(data)
    flux = subset_fluctuating_frequencies(compare_freqs)
    if len(flux) > 0:
        data[flux].plot(ax=axes[2], legend=False)
    # some fluctuating
    flux2 = subset_some_fluctuation(compare_freqs)
    for each in flux:
        flux2.remove(each)
    if len(flux2) > 0:
        data[flux2].plot(ax=axes[3], legend=False)
    # others
    oth = []
    for col in data.columns:
        if (col not in inc) & (col not in dec) & (col not in flux) & (col not in flux2):
            oth.append(col)
    if len(oth) > 0:
        data[oth].plot(ax=axes[4], legend=False)
    plt.savefig(plottitle + '_summary.png')
    plt.close()

def plot_line_individual(line, df, genome_id):
    data = df.loc[line.upper()].dropna(axis=1, how='all').fillna(0)
    inc, dec = subset_monotonic(data)
    compare_freqs = greater_than_last_gen(data)
    flux = subset_fluctuating_frequencies(compare_freqs)
    flux2 = subset_some_fluctuation(compare_freqs)
    for each in flux:
        flux2.remove(each)

    if len(inc) > 0:
        data[inc].plot(title = line.upper() + ':' + genome_id + ' is_monotonic_increasing', figsize=(20,15)).legend(loc='upper left', ncol=6)
        plt.savefig(line.upper() + '_' + genome_id + '_monotonic_increasing.png')
        plt.close()
    if len(dec) > 0:
        data[dec].plot(title = line.upper() + ':' + genome_id + ' is_monotonic_decreasing', figsize=(20,15)).legend(loc='upper right', ncol=6)
        plt.savefig(line.upper() + '_' + genome_id + '_monotonic_decreasing.png')
        plt.close()
    if len(flux) > 0:
        data[flux].plot(title = line.upper() + ':' + genome_id + ' fluctuating', figsize=(20,15)).legend(loc='upper right', ncol=6)
        plt.savefig(line.upper() + '_' + genome_id + '_fluctuating.png')
        plt.close()
    if len(flux2) > 0:
        data[flux2].plot(title = line.upper() + ':' + genome_id + ' some fluctuating', figsize=(20,15)).legend(loc='upper right', ncol=6)
        plt.savefig(line.upper() + '_' + genome_id + '_some_fluctuating.png')
        plt.close()

    oth = []
    for col in data.columns:
        if (col not in inc) & (col not in dec) & (col not in flux) & (col not in flux2):
            oth.append(col)
    if len(oth) > 0:
        data[oth].plot(title = line.upper() + ':' + genome_id + ' others', figsize=(20,15)).legend(loc='upper right', ncol=6)
        plt.savefig(line.upper() + '_' + genome_id + '_others.png')
        plt.close()
        
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


# frequencies from GD output
cp_gd = ComparePolymorphisms()
for line in evolution_lines:
    line_df = cp_gd.get_all_gd(line, input_directory, ancestor_gd_path)
    output_csv = output_directory + line + '_gd_freqs.csv'
    frequencies_df = cp_gd.gd_frequencies_to_df(line_df, save_csv=True, csv_filename=output_csv)
'''
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
    plot_line_individual(line, dvh_df, 'NC_002937')
    plot_line_summary(line, mm_df, (line + '-NC_005791'))
    plot_line_individual(line, mm_df, 'NC_005791')
    plot_line_summary(line, dvplasmid_df, (line + '-NC_005863'))
    plot_line_individual(line, dvplasmid_df, 'NC_005863')
    