import numpy as np
import pandas as pd
import tidy_breseq as tb

# Start data input
input_dir = '~/Repositories/dvh_mms2/depth_mle/'

d2_0_df = tb.get_mutations_df(input_dir + 'd2/annotated.gd', 'D2', '0', cov=True)
m1_0_df = tb.get_mutations_df(input_dir + 'm1/annotated.gd', 'M1', '0', cov=True)

df = pd.concat([d2_0_df, m1_0_df], ignore_index=True)
# End data input

df['name'] = df.culture.str.cat(df.generation.astype(str), sep='-')
names = df['name'].unique()

with open('samples.csv', 'w') as samples:
    samples.write('Name,Population,Generation\n')
    for name in names:
        population = name.split('-')[0]
        generation = name.split('-')[1]
        samples.write(name + ',' + population + ',' + generation + '\n')

with open('freqs.csv', 'w') as freqs:
    freqs.write('#CHROM,POS,REF,ALT,')
    count = 0
    while count < len(names) - 1:
        freqs.write(names[count] + ',')
        count += 1
    freqs.write(names[count] + '\n')
    
