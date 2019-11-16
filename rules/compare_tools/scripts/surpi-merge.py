import pandas as pd
import numpy as np
import logging

from ete3 import NCBITaxa

ncbi = NCBITaxa()

logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',datefmt='%d %b %Y %H:%M:%S',
                    filename=snakemake.log[0], level=logging.DEBUG)

virus_input=snakemake.input['vir']
bacteria_input=snakemake.input['bac']

vir_dic={}
for path in virus_input:
    sample=path.split('/')[1]
    vir_dic[sample]=path

bac_dic={}
for path in bacteria_input:
    sample=path.split('/')[1]
    bac_dic[sample]=path


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
    filter_counts = {k: v for k, v in dic.items() if v > threshold}
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


def names_to_taxid(ncbi,taxonomy, counts_dict):
    tax = {}
    target_ranks = ['superkingdom','phylum','order','family','genus','species']
    for i in taxonomy.keys():
        txid = list(ncbi.get_name_translator([taxonomy[i]]).values())
        if len(txid) == 1:
            taxid=txid[0][0]
            scientific_name = ncbi.translate_to_names([taxid])[0]
            tax[taxid] = {'taxid':str(taxid),'scientific_name': scientific_name, 'read_counts': counts_dict[i]}
            lineage = ncbi.get_lineage(taxid)
            names = ncbi.get_taxid_translator(lineage)
            ranks = ncbi.get_rank(lineage)
            rank_list = list(ranks.values())
            name_list = list(names.values())
            rank2names = dict(zip(rank_list, name_list))
            for sub_taxid in lineage:
                rank = ranks[sub_taxid]
                if rank in target_ranks:
                    tax[taxid][f"{rank}_taxid"] = sub_taxid
            for n, rank in enumerate(target_ranks):
                if rank in rank2names.keys():
                    tax[taxid][rank] = rank2names[rank]
                else:
                    previous_name = 'root'
                    for ranks_to_skip in range(1, len(target_ranks) + 1):
                        previous_rank = target_ranks[n - ranks_to_skip]
                        if previous_rank in rank2names.keys():
                            previous_name = rank2names[previous_rank]
                            break
                        if previous_rank not in rank2names.keys():
                            continue
                    tax[taxid][rank] = f'{previous_name}_{rank[0:1]}'
        if len(txid) > 1:
            logging.warning(f'More than one taxid matching input name: {taxonomy[i]} (at row number: {i})')
            tax = None
        if len(txid) == 0:
            logging.warning(f'No taxid matching for input name: {taxonomy[i]} (at row number: {i})')
            tax[i] = {'scientific_name': taxonomy[i], 'taxid': 'NA', 'read_counts': counts_dict[i]}
    return tax


def parse_surpi_table(table, ncbi, col_threshold, row_threshold):
    col_filtered_tb = filter_columns_low_read_counts(table, col_threshold)
    row_filtered_counts = get_filtered_row_counts(col_filtered_tb, row_threshold)
    filtered_row_taxonomy = get_taxonomy(col_filtered_tb, row_filtered_counts)
    dic = names_to_taxid(ncbi,filtered_row_taxonomy, row_filtered_counts)
    tb = pd.DataFrame.from_dict(dic, orient='index')
    tb = tb.replace(np.nan, 'NA')
    return tb


list_tables=[]

split_path=snakemake.output[0].split('/')
tool_path='/'.join(split_path[0:len(split_path)-1])

for sample in vir_dic.keys():
   v_tab=pd.read_csv(vir_dic[sample],sep='\t',low_memory=False)
   vir_f=parse_surpi_table(v_tab, ncbi, 5, 0)
   b_tab=pd.read_csv(bac_dic[sample],sep='\t',low_memory=False)
   bac_f=parse_surpi_table(b_tab, ncbi, 5, 0)
   full_tab = pd.concat([vir_f, bac_f],sort=False)
   full_tab = full_tab.replace(np.nan,'NA')
   full_tab=full_tab.groupby(['superkingdom','superkingdom_taxid','phylum','phylum_taxid','order','order_taxid','family','family_taxid','genus','genus_taxid','species','species_taxid','scientific_name','taxid']).sum()
   full_tab=full_tab.rename(columns={'read_counts':f'{sample}_counts'})
   full_tab[f'{sample}_percent'] = full_tab[f'{sample}_counts'].div(sum(full_tab[f'{sample}_counts'])).mul(100)
   full_tab.to_csv(f"{tool_path}/{sample}.tsv", sep='\t')
   list_tables.append(full_tab)

all_surpi_tab=pd.concat(list_tables,sort=False,axis=1)
all_surpi_tab.to_csv(snakemake.output[0],sep='\t')



