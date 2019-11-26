
import pandas

shortbred_outputs = snakemake.input 

all_samples = pandas.read_csv(snakemake.params["sample_table"], sep="\t", index_col=0)

o = open(snakemake.output[0], 'w')


for i, shortbred_output in enumerate(shortbred_outputs):
    sample_name = shortbred_output.split("/")[1]
    group_1 = all_samples.loc[sample_name, "group_1"]
    group_2 = all_samples.loc[sample_name, "group_2"]
    with open(shortbred_output, 'r') as f:
        for n, row in enumerate(f):
            if i == 0 and n == 0:
                o.write("sample_name\tgroup_1\tgroup_2\t" + row)
            if n != 0:
                o.write(f"{sample_name}\t{group_1}\t{group_2}\t" + row)
