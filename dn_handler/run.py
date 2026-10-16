import os

import pandas as pd

from dn_handler.input_handler import read_dn_df
from dn_handler.map_proteins import map_to_known_proteins

# CONFIG = {
#     'output_folder': 'output/te_dda',
#     'input_folder': 'test_data/transcription_errors',
#     'sample_info_file': 'test_data/transcription_errors/sample_info.csv',
#     'dn_method': 'casanovo',
#     'fasta_files': 'useful/UP000005640_9606.fasta',
#     'contam_fasta_files': 'useful/contam.fasta',
#     'score_cut_off': 0.8,
# }
CONFIG = {
    'output_folder': 'output/plasma_proteome',
    'input_folder': 'test_data/plasma_proteome',
    'sample_info_file': 'test_data/plasma_proteome/sample_info.csv',
    'dn_method': 'cascadia',
    'fasta_files': 'useful/UP000005640_9606.fasta',
    'contam_fasta_files': 'useful/contam.fasta',
    'score_cut_off': 0.8,
}


def main():
    os.makedirs(CONFIG['output_folder'], exist_ok=True)
    sample_info_df = pd.read_csv(CONFIG['sample_info_file'])

    dn_df, pivoted_df = read_dn_df(CONFIG['input_folder'], CONFIG['dn_method'], sample_info_df, CONFIG['score_cut_off'])
    dn_df.to_csv(f"{CONFIG['output_folder']}/complete_dn_df.csv", index=False)
    pivoted_df.to_csv(f"{CONFIG['output_folder']}/pivoted_dn_df.csv", index=False)

    pivoted_df = map_to_known_proteins(pivoted_df, CONFIG['fasta_files'], CONFIG['contam_fasta_files'])
    pivoted_df.to_csv(f"{CONFIG['output_folder']}/mapped_pivoted_dn_df.csv", index=False)
    pivoted_df = pd.read_csv(f"{CONFIG['output_folder']}/mapped_pivoted_dn_df.csv")

    unknown_df = pivoted_df[pivoted_df['unknown_mapping']].drop(columns=['contam_mapping','canonical_mapping','unknown_mapping'])


    return

if __name__ == '__main__':
    main()
