
CASANOVO_COLUMNS = [
    'sequence', 'search_engine_score[1]', 'spectra_ref',
]
CASANOVO_START_CODE = 'PSH\tsequence'
CASANOVO_RENAME_MAP = {
    'search_engine_score[1]': 'score',
    'spectra_ref': 'scan',
}

CASCADIA_COLUMNS = [
    'sequence', 'score', 'scan',
]