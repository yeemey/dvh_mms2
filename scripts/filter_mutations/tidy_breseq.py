import numpy as np
import pandas as pd
import re

def sum_cov(cov_str):
    cov_ints = cov_str.split('/')
    total_cov = 0
    for cov in cov_ints:
        total_cov += int(cov)
    return total_cov

def subset_gd_to_df(gdfile, sample_name, line, generation=np.nan, cov=False):
    '''
    Default returns one dataframe created from annotated.gd. All mutation rows are preserved,
    but only selected variables from each row, namely entry type, entry id,
    evidence id, genome id, position, mutation detail, frequency, and gene product.

    If cov=True, will return TWO dataframes, the first as above, the second reporting
    entry type, entry id, evidence id, genome id, position,reject reasons, prediction mode,
    polymorphism frequencies, major and minor coverages (i.e., major_cov, minor_cov),
    total coverage (total_cov), RA coverage (new_cov), JC coverage (new_junction_read_count),
    and MC-flanking coverage (left_outside_cov + right_outside_cov).
    '''
    df = pd.read_csv(gdfile, comment='#', names=range(200), dtype=str, sep='\t')
    df = df.dropna(axis=1, how='all')
    # https://stackoverflow.com/questions/27700591/reading-csv-files-with-messy-structure-with-pandas
    num_columns = len(df.columns)
    df.rename(columns = {0: 'entry_type', 1: 'entry_id', 2: 'evidence_id',
                         3: 'genome_id', 4: 'position', 5: 'mutation_detail'}, inplace=True)

    df_mutations = df[(df['entry_type'] == 'INS') | (df['entry_type'] == 'DEL') |
                      (df['entry_type'] == 'SNP') | (df['entry_type'] == 'SUB') |
                      (df['entry_type'] == 'MOB') | (df['entry_type'] == 'AMP') |
                      (df['entry_type'] == 'CON') | (df['entry_type'] == 'INV')].copy()

    for row in df_mutations.index:
        #check each column
        mut_col_index = 6
        while mut_col_index < num_columns:
            #1. mutation frequencies
            if re.match('frequency=', str(df_mutations.loc[row, mut_col_index])):
                df_mutations.loc[row, 'frequency'] = re.sub('frequency=', '', str(df_mutations.loc[row, mut_col_index]))
                if df_mutations.loc[row, 'frequency'] == 'NA':
                    df_mutations.loc[row, 'frequency'] = np.nan
            #2. gene products
            elif re.match('gene_product=', str(df_mutations.loc[row, mut_col_index])):
                df_mutations.loc[row, 'gene_product'] = re.sub('gene_product=', '', str(df_mutations.loc[row, mut_col_index]))
            mut_col_index += 1

    df_mutations = df_mutations[['entry_type', 'entry_id', 'evidence_id', 'genome_id',
                             'position', 'mutation_detail', 'frequency', 'gene_product']].copy()

    #insert sample name, line, and generation
    df_mutations.insert(0, 'sample', sample_name)
    df_mutations.insert(1, 'line', line)
    df_mutations.insert(2, 'generation', generation)
    #set frequencies type to float
    df_mutations['frequency'] = df_mutations['frequency'].astype(float)

    if cov == True:
        df_evidence = df[(df['entry_type'] == 'RA') | (df['entry_type'] == 'JC') |
                         (df['entry_type'] == 'MC') | (df['entry_type'] == 'UN')].copy()
        df_evidence.rename(columns = {6: 'REF', 7: 'ALT'}, inplace=True)

        for row in df_evidence.index:
            col_index = 8
            while col_index < num_columns:
                #3. polymorphism rejection reasons
                if re.match('reject=', str(df_evidence.loc[row, col_index])):
                    df_evidence.loc[row, 'reject'] = re.sub('reject=', '', str(df_evidence.loc[row, col_index]))
                #4. prediction type
                elif re.match('prediction=', str(df_evidence.loc[row, col_index])):
                    df_evidence.loc[row, 'prediction'] = re.sub('prediction=', '', str(df_evidence.loc[row, col_index]))
                #5. polymorphism mode frequencies
                elif re.match('polymorphism_frequency=', str(df_evidence.loc[row, col_index])):
                    df_evidence.loc[row, 'polymorphism_frequency'] = re.sub('polymorphism_frequency=', '', str(df_evidence.loc[row, col_index]))
                    if df_evidence.loc[row, 'polymorphism_frequency'] == 'NA':
                        df_evidence.loc[row, 'polymorphism_frequency'] = np.nan
                #6. major coverage counts
                elif re.match('major_cov=', str(df_evidence.loc[row, col_index])):
                    major_cov = re.sub('major_cov=', '', str(df_evidence.loc[row, col_index]))
                    df_evidence.loc[row, 'major_cov'] = sum_cov(major_cov)
                #7. minor coverage counts
                elif re.match('minor_cov', str(df_evidence.loc[row, col_index])):
                    minor_cov = re.sub('minor_cov=', '', str(df_evidence.loc[row, col_index]))
                    df_evidence.loc[row, 'minor_cov'] = sum_cov(minor_cov)
                #8. total coverage counts
                elif re.match('total_cov=', str(df_evidence.loc[row, col_index])):
                    total_cov = re.sub('total_cov=', '', str(df_evidence.loc[row, col_index]))
                    df_evidence.loc[row, 'total_cov'] = sum_cov(total_cov)
                #9. read alignment coverage counts
                elif re.match('new_cov=', str(df_evidence.loc[row, col_index])):
                    ra_cov = re.sub('new_cov=', '', str(df_evidence.loc[row, col_index]))
                    df_evidence.loc[row, 'ra_cov'] = sum_cov(ra_cov)
                #10. new junction coverage counts
                elif re.match('new_junction_read_count=', str(df_evidence.loc[row, col_index])):
                    df_evidence.loc[row, 'jc_cov'] = re.sub('new_junction_read_count=', '', str(df_evidence.loc[row, col_index]))
                #11. flanking coverage counts for missing coverage evidence
                elif re.match('left_outside_cov=', str(df_evidence.loc[row, col_index])):
                    left_cov = re.sub('left_outside_cov=', '', str(df_evidence.loc[row, col_index]))
                    if left_cov == 'NA':
                        left_cov = 0
                    else:
                        df_evidence.loc[row, 'left_cov'] = int(left_cov)
                elif re.match('right_outside_cov', str(df_evidence.loc[row, col_index])):
                    right_cov = re.sub('right_outside_cov=', '', str(df_evidence.loc[row, col_index]))
                    if right_cov == 'NA':
                        right_cov = 0
                    else:
                        df_evidence.loc[row, 'right_cov'] = int(right_cov)
                col_index += 1

        #set missing coverage col to 'NA' if no evidence
        if 'left_cov' in df_evidence.columns and 'right_cov' in df_evidence.columns:
            df_evidence[['left_cov', 'right_cov']].fillna(0)
            df_evidence['mc_cov'] = df_evidence.left_cov + df_evidence.right_cov
        else:
            df_evidence['mc_cov'] = np.nan
        #set reject col to 'NA' when no reject reason given.
        if 'reject' in df_evidence.columns:
            if (df_evidence.loc[row, 'reject'] == '') & (df_evidence.loc[row, 'evidence_id'] == '.'):
                df_evidence.loc[row, 'reject'] = np.nan
        else:
            df_evidence['reject'] = np.nan

        df_evidence = df_evidence[['entry_type', 'entry_id', 'genome_id', 'position', 'REF', 'ALT',
                     'reject', 'prediction', 'polymorphism_frequency', 'major_cov', 'minor_cov',
                     'total_cov', 'ra_cov', 'jc_cov', 'mc_cov']].copy()

        #insert sample name, line and generation
        df_evidence.insert(0, 'sample', sample_name)
        df_evidence.insert(1, 'line', line)
        df_evidence.insert(2, 'generation', generation)
        #set frequencies type to float
        df_evidence['polymorphism_frequency'] = df_evidence['polymorphism_frequency'].astype(float)

        return df_mutations, df_evidence

    else:
        return df_mutations

def combine_mutations_and_evidence(df_mutations, df_evidence):
    for evidence in df_mutations['evidence_id']:
        multi_evidence = evidence.split(',')
        count = 0
        while count < len(multi_evidence):
            mutation_row_index = df_mutations[df_mutations['evidence_id'] == evidence].index[0]
            df_mutations.loc[mutation_row_index, 'evidence_type'] = df_evidence.loc[int(multi_evidence[count])-1, 'entry_type']
            df_mutations.loc[mutation_row_index, 'REF'] = df_evidence.loc[int(multi_evidence[count])-1, 'REF']
            df_mutations.loc[mutation_row_index, 'ALT'] = df_evidence.loc[int(multi_evidence[count])-1, 'ALT']
            df_mutations.loc[mutation_row_index, 'reject'] = df_evidence.loc[int(multi_evidence[count])-1, 'reject']
            df_mutations.loc[mutation_row_index, 'prediction'] = df_evidence.loc[int(multi_evidence[count])-1, 'prediction']
            df_mutations.loc[mutation_row_index, 'polymorphism_frequency'] = df_evidence.loc[int(multi_evidence[count])-1, 'polymorphism_frequency']
            df_mutations.loc[mutation_row_index, 'major_cov'] = df_evidence.loc[int(multi_evidence[count])-1, 'major_cov']
            df_mutations.loc[mutation_row_index, 'minor_cov'] = df_evidence.loc[int(multi_evidence[count])-1, 'minor_cov']
            df_mutations.loc[mutation_row_index, 'total_cov'] = df_evidence.loc[int(multi_evidence[count])-1, 'total_cov']
            df_mutations.loc[mutation_row_index, 'ra_cov'] = df_evidence.loc[int(multi_evidence[count])-1, 'ra_cov']
            df_mutations.loc[mutation_row_index, 'jc_cov'] = df_evidence.loc[int(multi_evidence[count])-1, 'jc_cov']
            df_mutations.loc[mutation_row_index, 'mc_cov'] = df_evidence.loc[int(multi_evidence[count])-1, 'mc_cov']
            count += 1
    return df_mutations
