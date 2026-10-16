

def read_fasta(fasta_files):
    fasta_files = fasta_files.split(',')
    all_prots = {}
    for fasta_file in fasta_files:
        prot_id = None
        with open(fasta_file, 'r', ) as f: #TODO can be optimized
            for line in f:
                if line.startswith('>'):
                    prot_id = line[1:].strip().split()[0]
                    all_prots[prot_id] = ''
                else:
                    assert prot_id is not None, "FASTA format error: sequence line found before any header line"
                    all_prots[prot_id] += line.strip().replace('I', 'L')
    return all_prots

def map_to_known_proteins(pivoted_df, fasta_files, contam_fasta_files):
    # Implement the logic to map the de-novo dataframe to known proteins using the provided FASTA files
    all_prots = read_fasta(fasta_files)
    contam_prots = read_fasta(contam_fasta_files)
    contam_seqs = set(contam_prots.values())
    main_seqs = set(all_prots.values())

    pivoted_df['contam_mapping'] = pivoted_df['sequence'].apply(lambda seq: any(seq in s for s in contam_seqs))
    pivoted_df['canonical_mapping'] = pivoted_df['sequence'].apply(lambda seq: any(seq in s for s in main_seqs))
    pivoted_df['unknown_mapping'] = ~(pivoted_df['contam_mapping'] | pivoted_df['canonical_mapping'])

    return pivoted_df
    