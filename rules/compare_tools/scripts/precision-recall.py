import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

tables=snakemake.input.all_tools
rank=snakemake.params.rank
superkingdom=snakemake.params.superkingdom
gold_standard_file=snakemake.input.gold_standard
read_step=snakemake.params.read_step
max_iter=snakemake.params.max_read

gold_standard=pd.read_csv(gold_standard_file,sep='\t')
true_superkingdom=gold_standard[gold_standard['superkingdom']==superkingdom]
true_superkingdom.insert(loc=true_superkingdom.shape[1],column='tool',value=['gold_standard']*len(true_superkingdom))


tb_list=[]
tool_list=[]
for file in tables:
    tool = file.split('/')[1]
    tool_list.append(tool)
    tb=pd.read_csv(file,sep='\t')
    tb['tool']=[tool]*len(tb)
    tb_list.append(tb)

all_tools=pd.concat(tb_list,axis=0,sort=False)
tools_superkingdom=all_tools[all_tools.superkingdom==superkingdom]



def get_precision_recall_f1(true_tb, tool_tb, rank,tool_name,threshold):
    tool_tb=tool_tb[tool_tb['read_counts']>=threshold]
    true_list = list(set(true_tb[f'{rank}_taxid']))
    tool_list = list(set(tool_tb[f'{rank}_taxid']))
    tp=0
    fn=0
    for i in true_list:
        if tp<len(true_list) and i in true_list and i in tool_list:
            tp+=1
        if tp<len(true_list) and i in true_list and i not in tool_list:
            fn+=1
        if tp+fn==len(true_list):
            break

    fp = len(tool_list)-tp
    try:
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f1 = 2 * (precision * recall) / (precision + recall)
    except ZeroDivisionError:
        precision = 0
        recall = 0
        f1 = 0
    scores = {}
    scores[tool_name] = {'tool':tool_name,'tp': tp, 'fp': fp, 'fn': fn, 'precision': precision, 'recall': recall,
                            'F1_score': f1,'threshold':threshold}
    df = pd.DataFrame.from_dict(scores, orient='index')
    return df


def get_prc_table(true_tb, tool_tb, rank, tool_name, step, end):
    tables = []
    for threshold in range(0, end+step, step):
        tb = get_precision_recall_f1(true_tb, tool_tb, rank, tool_name, threshold)
        tables.append(tb)
    tables = pd.concat(tables)

    precisions = list(tables['precision'])
    recalls = list(tables['recall'])
    ap = [precisions[n] * (recalls[n] - recalls[n - 1]) for n in range(1, len(precisions))]#average precision
    aupr = round(abs(sum(ap)), 2)
    tables['AUPR'] = [aupr] * len(tables)
    return tables

gb_true=true_superkingdom.groupby(['tool','sample',f'{rank}_taxid'],as_index=False).sum()

tools_tables=[]
for tool in sorted(tool_list):
    if tool=='ezvir':
        gb_tools=tools_superkingdom.groupby(['tool','sample',f'{rank}_taxid'],as_index=False).sum()
        gb_tools['read_counts']=gb_tools['read_counts']/2
    else:
        gb_tools = tools_superkingdom.groupby(['tool', 'sample', f'{rank}_taxid'], as_index=False).sum()
        tools_tables.append(gb_tools)
gb_tools=pd.concat(tools_tables,sort=False)

max_threshold=max(gb_tools['read_counts'])#set the max range of precision and recall curves
if max_threshold>max_iter:#If the max number is too high, restrict it
    max_threshold=max_iter
subset=[]
for tool in sorted(tool_list):
        tool_tax=gb_tools[gb_tools['tool']==tool]
        table=get_prc_table(gb_true,tool_tax,rank,tool,read_step,max_threshold)
        subset.append(table)
subset=pd.concat(subset,sort=False)
sorted_subset=subset.sort_values(by=['AUPR','tool'],ascending=False)
sorted_subset.to_csv(snakemake.output.curve_table,sep='\t',index=False)

prc=sns.FacetGrid(data=sorted_subset,col='tool',col_wrap=3,hue='AUPR')
prc.map(plt.step,'recall','precision',where='pre',label='AUPR')
prc.map(plt.fill_between,'recall','precision',alpha=0.4,step='pre')
prc.add_legend()
prc.savefig(snakemake.output.curve)

indexed_subset=sorted_subset.set_index('threshold')
opt_read_thresh={}
for tool in tool_list:
    tool_tab=indexed_subset[indexed_subset.tool==tool]
    max_score=max(tool_tab.F1_score)
    index_list=list(tool_tab[tool_tab.F1_score==max_score].index)
    optimal_threshold=index_list[0]
    opt_read_thresh[tool]={'read_threshold':optimal_threshold,'max_F1_score':max_score}

opt_thresholds_tb=pd.DataFrame.from_dict(opt_read_thresh,orient='index')
opt_thresholds_tb.to_csv(snakemake.output.table,sep='\t')