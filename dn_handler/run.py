import os

import pandas as pd

from dn_handler.input_handler import read_dn_df
from dn_handler.map_proteins import map_to_known_proteins

# CONFIG = {
#     'output_folder': 'output/te_dda',
#     'input_folder': 'test_data/transcription_errors',
#     'sample_info_file': 'test_data/transcription_errors/sample_info.csv',
#     'dn_method': 'casanovo',
#     'fasta_files': 'resources/UP000005640_9606.fasta',
#     'contam_fasta_files': 'resources/contam.fasta',
#     'score_cut_off': 0.8,
# }
CONFIG = {
    'output_folder': 'output/plasma_proteome',
    'input_folder': 'test_data/plasma_proteome',
    'sample_info_file': 'test_data/plasma_proteome/sample_info.csv',
    'dn_method': 'cascadia',
    'fasta_files': 'resources/uniprot_swissprot_human.fasta',
    'contam_fasta_files': 'resources/crap.fasta',
    'score_cut_off': 0.8,
}
import subprocess

def add_trembl_isoforms(pivoted_df, output_folder):
    unknown_df = pivoted_df[pivoted_df['unknown_mapping']].drop(columns=['contam_mapping','canonical_mapping','unknown_mapping'])
    unknown_df.to_csv(f"{output_folder}/unknown_mapped_pivoted_dn_df.csv", index=False)
    unknown_df[['sequence']].to_csv(f"{output_folder}/unknown_sequences.csv", index=False, header=False)
    subprocess.run([
        'java', '-jar', 'resources/PeptideMatchCMD_1.1.jar', '-a', 'query',
        '-i', 'resources/uniprot_trembl_OX9606_index', '-Q',
        f"{CONFIG['output_folder']}/unknown_sequences.csv", '-l',
        '-e', '-o', f"{CONFIG['output_folder']}/unknown_sequences_mapped.txt"
    ])
    mapped_df = pd.read_csv(f"{CONFIG['output_folder']}/unknown_sequences_mapped.txt",
                            sep='\t', skiprows=2, header=None)
    mapped_df.columns = ['sequence', 'trembl_mapping']
    mapped_df = mapped_df.groupby('sequence', as_index=False).agg({'trembl_mapping': ';'.join})
    unknown_df = pd.merge(unknown_df, mapped_df, on='sequence', how='left')
    trembl_df = unknown_df[unknown_df['trembl_mapping'].notna()]
    unknown_df = unknown_df[unknown_df['trembl_mapping'].isna()].drop(columns=['trembl_mapping'])
    return trembl_df, unknown_df

def main():
    os.makedirs(CONFIG['output_folder'], exist_ok=True)
    sample_info_df = pd.read_csv(CONFIG['sample_info_file'])

    dn_df, pivoted_df = read_dn_df(CONFIG['input_folder'], CONFIG['dn_method'], sample_info_df, CONFIG['score_cut_off'])
    dn_df.to_csv(f"{CONFIG['output_folder']}/complete_dn_df.csv", index=False)
    pivoted_df.to_csv(f"{CONFIG['output_folder']}/pivoted_dn_df.csv", index=False)

    pivoted_df = map_to_known_proteins(pivoted_df, CONFIG['fasta_files'], CONFIG['contam_fasta_files'])

    trembl_df, unknown_df = add_trembl_isoforms(pivoted_df, CONFIG['output_folder'])

    unknown_df.to_csv(f"{CONFIG['output_folder']}/unknown_mapped_pivoted_dn_df.csv", index=False)
    trembl_df.to_csv(f"{CONFIG['output_folder']}/trembl_mapped_pivoted_dn_df.csv", index=False)


    return

if __name__ == '__main__':
    main()
