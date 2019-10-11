import pandas as pd
import numpy as np
from ete3 import NCBITaxa
ncbi = NCBITaxa()

virus_input=snakemake.input['vir']
bacteria_input=snakemake.input['bac']
sample_names=[name.split("/")[1] for name in virus_input]


def filter_columns_low_read_counts(table, threshold):
    colname_list = []
    for colname in table.columns:
        try:
            read_sum = int(sum(table[colname]))
            if read_sum >= threshold:
                colname_list.append(colname)
        except TypeError:  # If the column has strings
            colname_list.append(colname)
        except ValueError:  # If column has NaNs
            colname_list.append(colname)
    table_f = table[colname_list]
    table_f = table_f.replace('*', 'NA')
    table_f = table_f.replace(np.nan, 'NA')
    return table_f


def get_filtered_row_counts(filtered_col_tab, threshold):
    dic = dict(filtered_col_tab.sum(axis=1))
    filter_counts = {k: v for k, v in dic.items() if v >= threshold}
    return filter_counts


def get_taxonomy(filtered_col_tab, non_null_counts_rows):
    taxonomy = {}
    for i in non_null_counts_rows.keys():
        if filtered_col_tab.iloc[i]['species'] != 'NA':
            taxonomy[i] = filtered_col_tab.iloc[i]['species']
        elif filtered_col_tab.iloc[i]['species'] == 'NA' and filtered_col_tab.iloc[i]['genus'] != 'NA':
            taxonomy[i] = filtered_col_tab.iloc[i]['genus']
        elif filtered_col_tab.iloc[i]['species'] == 'NA' and filtered_col_tab.iloc[i]['genus'] == 'NA' and \
                filtered_col_tab.iloc[i]['family'] != 'NA':
            taxonomy[i] = filtered_col_tab.iloc[i]['family']
        else:
            taxonomy[i] = 'NA'
    return taxonomy


def names_to_taxid(taxonomy, counts_dict):
    names2taxid = {}
    target_ranks = ['species', 'genus', 'family', 'order', 'phylum', 'superkingdom']
    for i in taxonomy.keys():
        txid = list(ncbi.get_name_translator([taxonomy[i]]).values())
        if len(txid) == 1:
            lin_txid = ncbi.get_lineage(txid[0][0])
            lin_translation = ncbi.get_taxid_translator(lin_txid)
            taxid2rank = ncbi.get_rank(list(lin_translation.keys()))

            names2taxid[i] = {'name': taxonomy[i], 'taxid': txid[0][0], 'read_counts': counts_dict[i]}

            for k in taxid2rank.keys():
                if taxid2rank[k] in target_ranks:
                    names2taxid[i][taxid2rank[k]] = lin_translation[k]
        if len(txid) > 1:
            print(f'More than one taxid matching input name: {taxonomy[i]} (row index: {i})')
            names2taxid = None
        if len(txid) == 0:
            print(f'No taxid matching for input name: {taxonomy[i]} (row index: {i})')
            names2taxid[i] = {'name': taxonomy[i], 'taxid': 'NA', 'read_counts': counts_dict[i]}

    return names2taxid


def parse_surpi_table(table, col_threshold, row_threshold):
    col_filtered_tb = filter_columns_low_read_counts(table, col_threshold)
    row_filtered_counts = get_filtered_row_counts(col_filtered_tb, row_threshold)
    filtered_row_taxonomy = get_taxonomy(col_filtered_tb, row_filtered_counts)
    dic = names_to_taxid(filtered_row_taxonomy, row_filtered_counts)
    tb = pd.DataFrame.from_dict(dic, orient='index')
    tb = tb.replace(np.nan, 'NA')
    return tb

vir_dic=zip(dict(sample_names,virus_input))
bac_dic=zip(dict(sample_names,bacteria_input))
list_tables=[]
for i in sample_names:
   v_tab=pd.read_csv(vir_dic[i],sep='\t')
   vir_f=parse_surpi_table(v_tab,100,10)
   b_tab=pd.read_csv(bac_dic[i],sep='\t')
   bac_f=parse_surpi_table(b_tab,100,10)
   full_tab = pd.concat([vir_f, bac_f], sort=False)
   full_tab=full_tab.groupby(['taxid', 'name', 'species', 'genus', 'family', 'order', 'phylum', 'superkingdom']).sum()
   full_tab=full_tab.rename(columns={'read_counts':f'{i}'})
   list_tables.append(full_tab)

all_surpi_tab=pd.concat(list_tables,sort=False,axis=0,join='outer')
all_surpi_tab.to_csv(snakemake.output[0],sep='\t')



