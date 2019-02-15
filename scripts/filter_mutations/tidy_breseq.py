import numpy as np
import pandas as pd
import re

def sum_cov(cov_str):
    cov_ints = cov_str.split('/')
    total_cov = 0
    for cov in cov_ints:
        total_cov += int(cov)
    return total_cov

def subset_gd_to_df(gdfile, cov=False):
    '''
    Returns dataframe created from annotated.gd. All rows are preserved,
    but only selected variables from each row, namely entry type, entry id,
    evidence id, genome id, position, mutation detail, frequency, gene product,
    reject reasons, prediction mode, and polymorphism frequencies.

    If cov=True, will also report major and minor coverages (i.e., major_cov, minor_cov),
    total coverage (total_cov), RA coverage (new_cov), JC coverage (new_junction_read_count),
    and MC-flanking coverage (left_outside_cov + right_outside_cov).
    '''
    df = pd.read_csv(gdfile, comment='#', names=range(200), dtype=str, sep='\t')
    df = df.dropna(axis=1, how='all')
    # https://stackoverflow.com/questions/27700591/reading-csv-files-with-messy-structure-with-pandas
    num_columns = len(df.columns)
    df.rename(columns = {0: 'entry_type', 1: 'entry_id', 2: 'evidence_id',
                         3: 'genome_id', 4: 'position', 5: 'mutation_detail'}, inplace=True)
    for row in df.index:
        #check each column
        col_index = 6
        while col_index < num_columns:
            #1. mutation frequencies
            if re.match('frequency=', str(df.loc[row, col_index])):
                df.loc[row, 'frequency'] = re.sub('frequency=', '', str(df.loc[row, col_index]))
                if df.loc[row, 'frequency'] == 'NA':
                    df.loc[row, 'frequency'] = np.nan
            #2. gene products
            elif re.match('gene_product=', str(df.loc[row, col_index])):
                df.loc[row, 'gene_product'] = re.sub('gene_product=', '', str(df.loc[row, col_index]))
            #3. polymorphism rejection reasons
            elif re.match('reject=', str(df.loc[row, col_index])):
                df.loc[row, 'reject'] = re.sub('reject=', '', str(df.loc[row, col_index]))
            #4. prediction type
            elif re.match('prediction=', str(df.loc[row, col_index])):
                df.loc[row, 'prediction'] = re.sub('prediction=', '', str(df.loc[row, col_index]))
            #5. polymorphism mode frequencies
            elif re.match('polymorphism_frequency=', str(df.loc[row, col_index])):
                df.loc[row, 'polymorphism_frequency'] = re.sub('polymorphism_frequency=', '', str(df.loc[row, col_index]))
                if df.loc[row, 'polymorphism_frequency'] == 'NA':
                    df.loc[row, 'polymorphism_frequency'] = np.nan
            if cov == True:
                #6. major coverage counts
                if re.match('major_cov=', str(df.loc[row, col_index])):
                    major_cov = re.sub('major_cov=', '', str(df.loc[row, col_index]))
                    df.loc[row, 'major_cov'] = sum_cov(major_cov)
                #7. minor coverage counts
                elif re.match('minor_cov', str(df.loc[row, col_index])):
                    minor_cov = re.sub('minor_cov=', '', str(df.loc[row, col_index]))
                    df.loc[row, 'minor_cov'] = sum_cov(minor_cov)
                #8. total coverage counts
                elif re.match('total_cov=', str(df.loc[row, col_index])):
                    total_cov = re.sub('total_cov=', '', str(df.loc[row, col_index]))
                    df.loc[row, 'total_cov'] = sum_cov(total_cov)
                #9. read alignment coverage counts
                elif re.match('new_cov=', str(df.loc[row, col_index])):
                    ra_cov = re.sub('new_cov=', '', str(df.loc[row, col_index]))
                    df.loc[row, 'ra_cov'] = sum_cov(ra_cov)
                #10. new junction coverage counts
                elif re.match('new_junction_read_count=', str(df.loc[row, col_index])):
                    df.loc[row, 'jc_cov'] = re.sub('new_junction_read_count=', '', str(df.loc[row, col_index]))
                #11. flanking coverage counts for missing coverage evidence
                elif re.match('left_outside_cov=', str(df.loc[row, col_index])):
                    left_cov = re.sub('left_outside_cov=', '', str(df.loc[row, col_index]))
                    if left_cov == 'NA':
                        left_cov = 0
                    else:
                        df.loc[row, 'left_cov'] = int(left_cov)
                elif re.match('right_outside_cov', str(df.loc[row, col_index])):
                    right_cov = re.sub('right_outside_cov=', '', str(df.loc[row, col_index]))
                    if right_cov == 'NA':
                        right_cov = 0
                    else:
                        df.loc[row, 'right_cov'] = int(right_cov)
            col_index += 1
        #set reject col to 'NA' when no reject reason given.
        if 'reject' in df.columns:
            if (df.loc[row, 'reject'] == '') & (df.loc[row, 'evidence_id'] == '.'):
                df.loc[row, 'reject'] = np.nan
        else:
            df['reject'] = np.nan
    #set frequencies type to float
    df[['frequency', 'polymorphism_frequency']] = df[['frequency', 'polymorphism_frequency']].astype(float)
    if cov == True:
        df[['left_cov', 'right_cov']].fillna(0)
        df['mc_cov'] = df.left_cov + df.right_cov
        return df[['entry_type', 'entry_id', 'evidence_id', 'genome_id', 'position', 'mutation_detail',
                   'frequency', 'gene_product', 'reject', 'prediction', 'polymorphism_frequency',
                   'major_cov', 'minor_cov', 'total_cov', 'ra_cov', 'jc_cov', 'mc_cov']].copy()
    else:
        return df[['entry_type', 'entry_id', 'evidence_id', 'genome_id', 'position', 'mutation_detail',
                    'frequency', 'gene_product', 'reject', 'prediction', 'polymorphism_frequency']].copy()
    return df_subset

def select_mutation_rows(df):
    '''
    Selects mutation entry rows from dataframe.
    '''
    df_mutations = df[(df['entry_type'] == 'INS') |
                      (df['entry_type'] == 'DEL') |
                      (df['entry_type'] == 'SNP') |
                      (df['entry_type'] == 'SUB') |
                      (df['entry_type'] == 'MOB') |
                      (df['entry_type'] == 'AMP') |
                      (df['entry_type'] == 'CON') |
                      (df['entry_type'] == 'INV')]
    return df_mutations

def add_evidence_to_mutation_rows(df, cov=False):
    '''
    Adds new column for evidence.
    Evidence type taken from evidence entry rows, added to mutation entry rows.
    If cov=True, coverage counts are also added.
    '''
    df_mutations = select_mutation_rows(df)
    df['evidence'] = ''
    for evidence in df_mutations['evidence_id']:
        multi_evidence = evidence.split(',')
        count = 0
        while count < len(multi_evidence):
            entry_row_index = df_mutations[df_mutations['evidence_id'] == evidence].index[0]
            evidence_row_index = df[df['entry_id'] == multi_evidence[count]].index[0]
            df.loc[entry_row_index, 'prediction'] = df.loc[evidence_row_index, 'prediction']
            df.loc[entry_row_index, 'polymorphism_frequency'] = df.loc[evidence_row_index, 'polymorphism_frequency']
            df.loc[entry_row_index, 'evidence'] = df.loc[entry_row_index, 'evidence'] + df.loc[evidence_row_index, 'entry_type']
            if cov == True:
                df.loc[entry_row_index, 'major_cov'] = df.loc[evidence_row_index, 'major_cov']
                df.loc[entry_row_index, 'minor_cov'] = df.loc[evidence_row_index, 'minor_cov']
                df.loc[entry_row_index, 'total_cov'] = df.loc[evidence_row_index, 'total_cov']
                df.loc[entry_row_index, 'ra_cov'] = df.loc[evidence_row_index, 'ra_cov']
                df.loc[entry_row_index, 'jc_cov'] = df.loc[evidence_row_index, 'jc_cov']
                df.loc[entry_row_index, 'mc_cov'] = df.loc[evidence_row_index, 'mc_cov']
            count += 1
    return df

def get_mutations_df(gdfile, culture_name, generation=np.nan, cov=False):
    '''
    Adds new columns for culture name, and sampled generation.
    Returns df of only mutation entry rows, with added evidence column.
    If cov=True, coverage counts are also added.
    '''
    if cov == True:
        gd_subset = subset_gd_to_df(gdfile, cov=True)
        gd_subset_df = add_evidence_to_mutation_rows(gd_subset, cov=True)
    else:
        gd_subset = subset_gd_to_df(gdfile)
        gd_subset_df = add_evidence_to_mutation_rows(gd_subset)
    gd_subset_df.insert(0, 'culture', culture_name)
    gd_subset_df.insert(1, 'generation', generation)
    mutations_df = select_mutation_rows(gd_subset_df).copy()
    return mutations_df
