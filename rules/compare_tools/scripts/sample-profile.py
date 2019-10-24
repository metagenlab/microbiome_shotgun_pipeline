from datetime import date
import pandas as pd
today=date.today()
d=today.strftime("%d/%m/%Y")
sample=snakemake.wildcards.sample
table=pd.read_csv(snakemake.input[0],sep='\t')
table=table[['last_taxid','last_rank','taxid_path',f'{sample}_percent']]
table.columns=['TAXID','RANK','TAXPATH','PERCENTAGE']
str_tb=table.to_string(index=False)
profile_file=open(snakemake.output[0],"w+")
profile_file.write(f'@SampleID:{sample}\n\
@Version:0.9.1\n\
@Ranks:superkingdom|phylum|order|family|genus|species\n\
@TaxonomyID:ncbi-taxonomy_{d}\n\
@@{str_tb}')
profile_file.close()


# def open_profile_from_tsv(file_path, normalize):
#     header = {}
#     column_name_to_index = {}
#     profile = []
#     samples_list = []
#     predictions_dict = {}
#     reading_data = False
#     got_column_indices = False
#
#     with open(file_path) as read_handler:
#         for line in read_handler:
#             if len(line.strip()) == 0 or line.startswith("#"):
#                 continue
#             line = line.rstrip('\n')
#
#             # parse header with column indices
#             if line.startswith("@@"):
#                 for index, column_name in enumerate(line[2:].split('\t')):
#                     column_name_to_index[column_name] = index
#                 index_rank, index_taxid, index_percentage, index_taxpath, index_taxpathsn = get_column_indices(column_name_to_index)
#                 got_column_indices = True
#                 reading_data = False
#                 continue
#
#             # parse header with metadata
#             if line.startswith("@"):
#                 # if last line contained sample data and new header starts, store profile for sample
#                 if reading_data:
#                     if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
#                         if len(profile) > 0:
#                             samples_list.append((header['SAMPLEID'], header, profile))
#                             profile = []
#                             predictions_dict = {}
#                     else:
#                         logging.getLogger('opal').critical("Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.\n".format(file_path))
#                         raise RuntimeError
#                     header = {}
#                 reading_data = False
#                 got_column_indices = False
#                 key, value = line[1:].split(':', 1)
#                 header[key.upper()] = value.strip()
#                 continue
#
#             if not got_column_indices:
#                 logging.getLogger('opal').critical("Header line starting with @@ in file {} is missing or at wrong position.\n".format(file_path))
#                 raise RuntimeError
#
#             reading_data = True
#             row_data = line.split('\t')
#
#             taxid = row_data[index_taxid]
#             # if there is already a prediction for taxon, only sum abundance
#             if taxid in predictions_dict:
#                 prediction = predictions_dict[taxid]
#                 prediction.percentage += float(row_data[index_percentage])
#             else:
#                 prediction = Prediction()
#                 predictions_dict[taxid] = prediction
#                 prediction.taxid = row_data[index_taxid]
#                 prediction.rank = row_data[index_rank]
#                 prediction.percentage = float(row_data[index_percentage])
#                 prediction.taxpath = row_data[index_taxpath]
#                 if isinstance(index_taxpathsn, int):
#                     prediction.taxpathsn = row_data[index_taxpathsn]
#                 else:
#                     prediction.taxpathsn = None
#                 profile.append(prediction)
#
#     # store profile for last sample
#     if 'SAMPLEID' in header and 'VERSION' in header and 'RANKS' in header:
#         if reading_data and len(profile) > 0:
#             samples_list.append((header['SAMPLEID'], header, profile))
#     else:
#         logging.getLogger('opal').critical("Header in file {} is incomplete. Check if the header of each sample contains at least SAMPLEID, VERSION, and RANKS.\n".format(file_path))
#         raise RuntimeError
#
#     if normalize:
#         normalize_samples(samples_list)
#
#     return samples_list