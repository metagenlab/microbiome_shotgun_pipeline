import pandas as pd
from Bio import Entrez
Entrez.email = 'farid.chaabane@unil.ch'
Entrez.api_key ='4c437af76fccca44ae0f28ba232ae6c91909'
tab_ganon=pd.read_csv(snakemake.input[0],sep='\t',header=None,dtype='str')
tab_ganon.columns=['SeqID','length','taxid','AssemblyID']
na_acc_list=tab_ganon[tab_ganon.length=='na']['AssemblyID']#Extract assembly accessions for rows with na sequence length
na_list_index=list(na_acc_list.index)
na_list_index_dic=dict(zip(na_acc_list,na_list_index))
len_dic={}
for i in na_acc_list:
    uid=Entrez.read(Entrez.esearch(db='assembly', term=i))['IdList']
    assembly_info=Entrez.read(Entrez.esummary(db='assembly', id=uid))
    gen_len=assembly_info['DocumentSummarySet']['DocumentSummary'][0]['Meta'].split('Stat category="total_length" sequence_tag="all">')[1].rsplit('</Stat>')[0]#Extract sequence len from NCBI metadata
    len_dic[i]=gen_len
    tab_ganon.iloc[na_list_index_dic[i]]['length'] = len_dic[i]

tab_ganon.to_csv(snakemake.output[0],sep='\t',header=None,index=None)