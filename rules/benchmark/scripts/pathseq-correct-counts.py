from ete3 import NCBITaxa
import pandas as pd
import copy
ncbi = NCBITaxa()


path=snakemake.input[0]
taxid2count= pd.read_csv(path,delimiter='\t').set_index("tax_id").to_dict()["unambiguous"]
taxid2name= pd.read_csv(path,delimiter='\t').set_index("tax_id").to_dict()["name"]
taxid_set = set(taxid2count.keys())


# check if possible to retrieve species name
# if not, discard taxid
keep = []
for taxid in taxid_set:
    conv = ncbi.translate_to_names([taxid])[0]
    if str(conv) == str(taxid) and taxid2count[taxid]==0:
        print(f"problem with {taxid}: {taxid2name[taxid]}")
        print(f'read counts: {taxid2count[taxid]}')
        print('Discarding entry')
    elif str(conv) == str(taxid) and taxid2count[taxid]>0:
        print(f"problem with {taxid}: {taxid2name[taxid]}")
        print(f'read counts: {taxid2count[taxid]}')
        print('Discarding entry')
    else:
        keep.append(taxid)

taxid2corrected_counts = {}
tree = ncbi.get_topology(keep)

# ATTENTION: ete add taxid absent from the initial dataset ==> try to check if in dico
n = 0
for node in tree.traverse():
    n+=1
    counts = []
    for i in node.children:
        try:
            counts.append(taxid2count[int(i.name)])
        except KeyError:
            #print("problem with %s" % node.name)
            pass
    try:
        taxid2corrected_counts[node.name] = taxid2count[int(node.name)] - sum(counts)
    except:
        #print("problem with %s" % node.name)
        pass

taxid_info=pd.read_csv(path,sep='\t').set_index('tax_id').to_dict(orient='index')
new_scores=copy.deepcopy(taxid_info)

for taxid in new_scores.keys():
    if taxid!=1:
        try:
            new_scores[taxid]['unambiguous']=taxid2corrected_counts[f'{taxid}']
        except KeyError:
            new_scores[taxid]['unambiguous']=0#Discard entries for which we cannot correct cumulated counts
    elif taxid==1:
        new_scores[taxid]['unambiguous'] = 0

tb=pd.DataFrame.from_dict(new_scores,orient='index')
tb.to_csv(snakemake.output[0],sep='\t')