from sys import argv
import numpy as np
import pandas as pd
import tidy_breseq as tb

script, line, input_dir, sample, generation = argv
#e.g., line=UA3, sample=U1-15, generation=100

dfmut, dfev = tb.subset_gd_to_df(input_dir + '/output/evidence/annotated.gd', sample, line, generation, cov=True)
complete_df = tb.combine_mutations_and_evidence(dfmut, dfev)
complete_df.to_csv(sample+'_mutev.csv', index=False)
dfmut.to_csv(sample+'_mutations.csv', index=False)
dfev.to_csv(sample+'_evidence.csv', index=False)
