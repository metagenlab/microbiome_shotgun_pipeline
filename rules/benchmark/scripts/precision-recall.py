import pandas as pd

###Snakemake inputs, wildcards and params
file=snakemake.input.tool
gold_standard_file=snakemake.input.gold_standard
rank=snakemake.params.rank
superkingdom=snakemake.params.superkingdom
read_step=snakemake.params.read_step
max_iter=snakemake.params.max_read
tool = snakemake.wildcards.tool

#Gold standard table
gold_standard=pd.read_csv(gold_standard_file,sep='\t')
true_superkingdom=gold_standard[gold_standard['superkingdom']==superkingdom]
true_superkingdom.insert(loc=true_superkingdom.shape[1],column='tool',value=['gold_standard']*len(true_superkingdom))
true_tb=true_superkingdom.groupby(['tool','sample',f'{rank}_taxid'],as_index=False).sum()


# Tool table
tool_tb=pd.read_csv(file,sep='\t')
tool_tb=tool_tb[tool_tb.superkingdom==superkingdom]
tool_tb['tool']=[tool]*len(tool_tb)
tool_tb=tool_tb.groupby(['tool', 'sample', f'{rank}_taxid'], as_index=False).sum()
if tool=='ezvir' or tool=='metaphlan':
    tool_tb['read_counts']=tool_tb['read_counts']//2


def get_precision_recall_f1_per_sample(true_tb, tool_tb, rank, tool_name, threshold):
    target_col = f'{rank}_taxid'
    if rank == 'scientific_name':
        target_col = 'taxid'
    tool_tb = tool_tb[tool_tb['read_counts'] >= threshold]

    samples = []
    for sample in list(set(true_tb['sample'])):
        true_list = list(set(true_tb[true_tb['sample'] == sample][target_col]))
        true_list = [str(int(float(taxid))) for taxid in true_list]  # for converting float taxid to integers

        tool_list = list(set(tool_tb[tool_tb['sample'] == sample][target_col]))
        tool_list = [str(int(float(taxid))) for taxid in tool_list]  # for converting float taxid to integers
        tp = 0
        fn = 0
        for i in true_list:
            if tp < len(true_list) and i in true_list and i in tool_list:
                tp += 1
            if tp < len(true_list) and i in true_list and i not in tool_list:
                fn += 1
            if tp + fn == len(true_list):
                break

        fp = len(tool_list) - tp
        try:
            precision = tp / (tp + fp)
            recall = tp / (tp + fn)
            f1 = 2 * (precision * recall) / (precision + recall)
        except ZeroDivisionError:
            precision = 0
            recall = 0
            f1 = 0
        scores = {}
        scores[tool_name] = {'tool': tool_name, 'sample': sample, 'tp': tp, 'fp': fp, 'fn': fn, 'precision': precision,
                             'recall': recall,
                             'F1_score': f1, 'threshold': threshold}
        df = pd.DataFrame.from_dict(scores, orient='index')
        samples.append(df)
    tab = pd.concat(samples)
    return tab


def get_prc_table(true_tb, tool_tb, rank, tool_name, step, end):
    tables = []
    for threshold in range(0, end+step, step):
        tb = get_precision_recall_f1_per_sample(true_tb, tool_tb, rank, tool_name, threshold)
        tables.append(tb)
    tables = pd.concat(tables)
    tables=tables.sort_values(by=['sample','threshold'])
    samptabs=[]
    for sample in sorted(list(set(tables['sample']))):
        samptab=tables[tables['sample']==sample]
        precisions = list(samptab['precision'])
        recalls = list(samptab['recall'])
        ap = [precisions[n] * (recalls[n] - recalls[n - 1]) for n in range(1, len(precisions))]#average precision
        aupr = round(abs(sum(ap)), 2)
        samptab.insert(loc=samptab.shape[1],column='AUPR',value=[aupr] * len(samptab))
        samptabs.append(samptab)
    auprtab=pd.concat(samptabs)
    return auprtab

max_threshold=max(tool_tb['read_counts'])#set the max range of precision and recall curves
if max_threshold>max_iter:#If the max number is too high, restrict it
    max_threshold=max_iter
aupr_tab=get_prc_table(true_tb,tool_tb,rank,tool,read_step,max_threshold)
aupr_tab.to_csv(snakemake.output[0],sep='\t',index=False)








