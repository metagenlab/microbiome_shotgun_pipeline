import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


tables=snakemake.input.all_tools
rank=snakemake.params.rank
superkingdom=snakemake.params.superkingdom
gold_standard=snakemake.input.gold_standard

true_superkingdom=gold_standard[gold_standard['superkingdom']==superkingdom]
max_threshold=max(true_superkingdom.read_counts)
read_step=snakemake.params.read_step

tb_list=[]
for file in tables:
    tool = file.split('/')[1]
    tb=pd.read_csv(file,sep='\t')
    tb['tool']=[tool]*len(tb)
    tb_list.append(tb)

all_tools=pd.concat(tb_list,axis=0,sort=False)
tools_superkingdom=all_tools[all_tools.superkingdom==superkingdom]



def get_precision_recall_f1(true_tb, tool_tb, rank, tool_name):
    if rank=='scientific_name':
        taxid_col='taxid'
    else:
        taxid_col=f'{rank}_taxid'
    true_list = list(true_tb.groupby([taxid_col], as_index=False).sum()[taxid_col])
    tool_list = list(tool_tb.groupby([taxid_col], as_index=False).sum()[taxid_col])
    tp_list = [i for i in true_list if i in true_list and i in tool_list]
    fp_list = [i for i in tool_list if i not in true_list and i in tool_list]
    fn_list = [i for i in true_list if i in true_list and i not in tool_list]
    tp = len(tp_list)
    fp = len(fp_list)
    fn = len(fn_list)
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
                            'F1_score': f1}
    df = pd.DataFrame.from_dict(scores, orient='index')
    return df


def get_threshold_table(true_tb, all_tools, max_read, step, rank):
    tool_list = list(set(all_tools['tool']))
    thresholds = list(range(0, max_read, step))
    tables = []
    for threshold in thresholds:
        for tool_name in tool_list:
            tool_tb = all_tools[all_tools['tool'] == tool_name]
            tb = get_precision_recall_f1(true_tb, tool_tb[tool_tb['read_counts'] >= threshold], rank, tool_name)
            tb['threshold'] = threshold
            tables.append(tb)
    pr_tb = pd.concat(tables)
    return pr_tb

def get_auprc(table):
    precision=list(table['precision'])
    recall=list(table['recall'])
    ap=[]
    for n in range(1,len(precision)):
        ap.append(precision[n]*(recall[n]-recall[n-1]))
    auc=round(abs(sum(ap)),2)
    return auc

pr_table=get_threshold_table(true_superkingdom,tools_superkingdom,max_threshold+read_step,read_step,rank)

pos_val_tb=pr_table[pr_table.F1_score!=0]

subset=[]
tool_list=list(set(pos_val_tb['tool']))
for tool in tool_list:
        table=pos_val_tb[pos_val_tb['tool']==tool]
        auc=get_auprc(table)
        table['AUC']=[auc]*len(table)
        subset.append(table)
subset=pd.concat(subset,sort=False)
sorted_subset=subset.sort_values(by=['AUC','tool'],ascending=False)

roc=sns.FacetGrid(data=sorted_subset,col='tool',col_wrap=4,hue='AUC')
roc.map(plt.step,'recall','precision',where='pre',label='AUC')
roc.map(plt.fill_between,'recall','precision',alpha=0.4,step='pre')
roc.add_legend()
roc.savefig(snakemake.output.curve)

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