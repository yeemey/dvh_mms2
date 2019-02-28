import numpy as np
import pandas as pd
import tidy_breseq as tb

# Start data input
input_dir = '~/Repositories/dvh_mms2/depth_mle/'
d2mut, d2ev = tb.subset_gd_to_df(input_dir + 'd2/annotated.gd', 'D2-0', 'D2', '0', cov=True)
m1mut, m1ev = tb.subset_gd_to_df(input_dir + 'm1/annotated.gd', 'M1-0', 'M1', '0', cov=True)

d2 = tb.combine_mutations_and_evidence(d2mut, d2ev)
m1 = tb.combine_mutations_and_evidence(m1mut, m1ev)

df = pd.concat([d2, m1], ignore_index=True)
# End data input

with open('samples.csv', 'w') as samples:
    samples.write('Name,Population,Generation\n')
    samples_df = df[['sample', 'line', 'generation']].drop_duplicates()
    for row in samples_df.index:
        samples.write(samples_df.loc[row, 'sample'] + ',' + samples_df.loc[row, 'line'] + ',' + str(samples_df.loc[row, 'generation'] + '\n'))

df_RA = df[df['evidence_type'] == 'RA']

freqs_df = df_RA[['genome_id', 'position', 'REF', 'ALT', 'sample', 'polymorphism_frequency']]
pvt_freqs = freqs_df.pivot_table(index=['genome_id', 'position', 'REF', 'ALT'], columns='sample', values='polymorphism_frequency')
with open('freqs.csv', 'w') as freqs:
    header = '#CHROM,POS,REF,ALT,'
    count = 0
    while count < len(pvt_freqs.columns) - 1:
        header += pvt_freqs.columns[count] + ','
        count += 1
    header += pvt_freqs.columns[count] + '\n'
    freqs.write(header)
    for row in pvt_freqs.index:
        freqs.write(row[0] + ',' + row[1] + ',' + row[2] + ',' + row[3] + ',')
        sample_count = 0
        while sample_count < len(pvt_freqs.columns) - 1:
            freqs.write(str(pvt_freqs.loc[row, pvt_freqs.columns[sample_count]]) + ',')
            sample_count += 1
        freqs.write(str(pvt_freqs.loc[row, pvt_freqs.columns[sample_count]]) + '\n')

cov_df = df_RA[['genome_id', 'position', 'REF', 'ALT', 'sample', 'ra_cov']]
pvt_cov = cov_df.pivot_table(index=['genome_id', 'position', 'REF', 'ALT'], columns='sample', values='ra_cov').fillna(0)
with open('coverage.csv', 'w') as cov:
    header = '#CHROM,POS,REF,ALT,'
    count = 0
    while count < len(pvt_cov.columns) - 1:
        header += pvt_cov.columns[count] + ','
        count += 1
    header += pvt_cov.columns[count] + '\n'
    cov.write(header)
    for row in pvt_cov.index:
        cov.write(row[0] + ',' + row[1] + ',' + row[2] + ',' + row[3] + ',')
        sample_count = 0
        while sample_count < len(pvt_cov.columns) - 1:
            cov.write(str(pvt_cov.loc[row, pvt_cov.columns[sample_count]]) + ',')
            sample_count += 1
        cov.write(str(pvt_cov.loc[row, pvt_cov.columns[sample_count]]) + '\n')
