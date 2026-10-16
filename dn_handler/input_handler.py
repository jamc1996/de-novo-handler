
import os

import pandas as pd
from dn_handler.constants import (
    CASANOVO_COLUMNS, CASANOVO_START_CODE, CASANOVO_RENAME_MAP,
    CASCADIA_COLUMNS
)


def read_dn_df(input_folder, dn_method, sample_info_df, score_cut_off):
    if dn_method == 'casanovo':
        dn_df = read_casanovo_dfs(input_folder, sample_info_df, score_cut_off)
    elif dn_method == 'cascadia':
        dn_df = read_cascadia_dfs(input_folder, sample_info_df, score_cut_off)

    dn_df = dn_df.groupby(
        ['sample_name', 'sequence'], as_index=False
    ).agg({'score': 'max', 'scan': list})
    dn_df = dn_df.rename(columns={'score': 'max_score', 'scan': 'scans'})
    dn_df['scans'] = dn_df['scans'].apply(lambda x: ' '.join(sorted(x)))

    dn_df['present'] = 1
    pivoted_df = dn_df.pivot_table(
        index='sequence',
        columns='sample_name',
        values='present',
        aggfunc='max',
        fill_value=0,
    ).astype(int)
    pivoted_df = pivoted_df.reset_index().rename(columns={'index': 'sequence'})

    return dn_df,pivoted_df

def clean_up_columns(casa_df, sample_info, score_cut_off):

    casa_df = casa_df[casa_df['score'] >= score_cut_off]
    casa_df['sequence'] = casa_df['sequence'].apply(lambda x: x.replace('I', 'L'))

    casa_df['probe_name'] = sample_info['probe_name']
    casa_df['sample_name'] = sample_info['sample_name']
    casa_df['scan'] = casa_df['scan'].apply(lambda x : f"{sample_info['probe_name']}_{x}")

    # remove all characters that are not upper case letters
    casa_df['sequence'] = casa_df['sequence'].apply(lambda x: ''.join([c for c in x if c.isupper()]))
    return casa_df

def read_casanovo_dfs(input_folder, sample_info_df, score_cut_off):
    sample_info_df['probe_path'] = sample_info_df['raw_probe_name'].apply(
        lambda x: f"{input_folder}/{x}.mztab"
    )

    all_casa_dfs = []
    for _, df_row in sample_info_df.iterrows():
        skip_rows = 0
        with open(df_row['probe_path'], 'r') as f:
            for line in f:
                if line.startswith(CASANOVO_START_CODE):
                    break
                skip_rows += 1
        casa_df = pd.read_csv(
            df_row['probe_path'], sep='\t', skiprows=skip_rows,
            usecols=CASANOVO_COLUMNS
        )
        casa_df = casa_df.rename(columns=CASANOVO_RENAME_MAP)
        casa_df = clean_up_columns(casa_df, df_row, score_cut_off)
        all_casa_dfs.append(casa_df)

    return pd.concat(all_casa_dfs, ignore_index=True)
        

def read_cascadia_dfs(input_folder, sample_info_df, score_cut_off):
    sample_info_df['probe_path'] = sample_info_df['raw_probe_name'].apply(
        lambda x: f"{input_folder}/{x}.ssl"
    )

    all_cascadia_dfs = []
    for _, df_row in sample_info_df.iterrows():
        cascadia_df = pd.read_csv(
            df_row['probe_path'], sep='\t', usecols=CASCADIA_COLUMNS
        )
        cascadia_df = clean_up_columns(cascadia_df, df_row, score_cut_off)

        all_cascadia_dfs.append(cascadia_df)

    return pd.concat(all_cascadia_dfs, ignore_index=True)
