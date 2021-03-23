import pandas as pd
import numpy as np

readlen = snakemake.params.readlen
rank = snakemake.params.rank
pairing = snakemake.params.pairing
taxtb = pd.read_csv(snakemake.input.tab, sep='\t')
prok = pd.read_csv(snakemake.input.prok, sep='\t', low_memory=False)
euk = pd.read_csv(snakemake.input.euk, sep='\t')
vir = pd.read_csv(snakemake.input.vir, sep='\t')
tables = [euk, vir, prok]


def get_sizes(tb):
    """
    function that returns taxids and their mean genome sizes
    """
    if 'Size (Mb)' in tb.columns:
        sizecol = 'Size (Mb)'
        factor = 10**6
    elif 'Size (Kb)' in tb.columns:
        sizecol = 'Size (Kb)'
        factor = 10**3
    tb = tb.set_index('Status').drop(['Contig', 'Scaffold'])  # Drop Contigs and scafflods
    tb.replace('-', np.nan, inplace=True)
    tb = tb.astype({sizecol: 'float64'})
    tb.reset_index(inplace=True)
    subset = tb[['TaxID', sizecol]].groupby('TaxID', as_index=False).mean()
    subset.insert(loc=subset.shape[1], column='meansize', value=subset[sizecol].mul(factor))
    subset.drop(columns=sizecol, inplace=True)
    return subset


cat = []
for t in tables:
    cat.append(get_sizes(t))
all_sizes = pd.concat(cat)  # Table with a list of taxids with their mean genome size
all_sizes = all_sizes.astype({'TaxID': 'int64'})
s = taxtb[(taxtb[f'{rank}_taxid'] != 'na')]  # select only hits identified at the species level for example
gb = s.groupby(['superkingdom', f'{rank}', f'{rank}_taxid', 'sample'], as_index=False)['read_counts'].sum()
gb = gb.rename(columns={f'{rank}_taxid': 'TaxID'})
gb = gb.astype({'TaxID': 'int64'})
merged = gb.merge(all_sizes, on='TaxID', how='left').replace(np.nan, 0)
merged.insert(loc=merged.shape[1], column='coverage',
              value=round(merged['read_counts']*readlen*pairing/merged['meansize']))
merged.replace(np.inf, 0, inplace=True)
merged.sort_values(by='coverage', ascending=False).to_csv(snakemake.output[0], sep='\t', index=None)
