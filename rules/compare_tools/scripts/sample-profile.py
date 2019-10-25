from datetime import date
import pandas as pd
today=date.today()
d=today.strftime("%d.%m.%Y")
sample=snakemake.wildcards.sample
tool=snakemake.wildcards.tool
table=pd.read_csv(snakemake.input[0],sep='\t')
table=table[['last_taxid','last_rank','taxid_path',f'{sample}_percent']]
table.columns=['TAXID','RANK','TAXPATH','PERCENTAGE']

row_list=[]
for index, row in table.iterrows():
    taxid=row['TAXID']
    rank=row['RANK']
    path=row['TAXPATH']
    percent=float(row['PERCENTAGE'])
    row_list.append(f"{taxid}\t{rank}\t{path}\t{percent}")

str_tb='\n'.join(row_list)

profile_file=open(snakemake.output[0],"w+")
profile_file.write(f'@SampleID:{sample}\n\
@Version:0.9.1\n\
@Ranks:superkingdom|phylum|order|family|genus|species\n\
@TaxonomyID:ncbi-taxonomy_{d}\n\
@__program__:{tool}\n\
@@TAXID\tRANK\tTAXPATH\tPERCENTAGE\n\
{str_tb}')
profile_file.close()


