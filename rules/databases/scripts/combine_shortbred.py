
import pandas

shortbred_outputs = snakemake.input["shortbred_output"]
silix_98 = snakemake.input["silix_98"]
silix_95 = snakemake.input["silix_95"]
silix_90 = snakemake.input["silix_90"]
silix_80 = snakemake.input["silix_80"]

accession2silix_98 = pandas.read_csv(silix_98,
                                     delimiter='\t',
                                     names=["silix_name","accession"]).set_index("accession").to_dict()["silix_name"]

accession2silix_95 = pandas.read_csv(silix_95,
                                     delimiter='\t',
                                     names=["silix_name","accession"]).set_index("accession").to_dict()["silix_name"]

accession2silix_90 = pandas.read_csv(silix_90,
                                     delimiter='\t',
                                     names=["silix_name","accession"]).set_index("accession").to_dict()["silix_name"]

accession2silix_80 = pandas.read_csv(silix_80,
                                     delimiter='\t',
                                     names=["silix_name","accession"]).set_index("accession").to_dict()["silix_name"]

all_samples = pandas.read_csv(snakemake.params["sample_table"], sep="\t", index_col=0)

o = open(snakemake.output[0], 'w')


for i, shortbred_output in enumerate(shortbred_outputs):
    sample_name = shortbred_output.split("/")[1]
    group_1 = all_samples.loc[sample_name, "group_1"]
    group_2 = all_samples.loc[sample_name, "group_2"]
    with open(shortbred_output, 'r') as f:
        for n, row in enumerate(f):
            if i == 0 and n == 0:
                o.write("sample_name\tgroup_1\tgroup_2\silix_98\silix_95\silix_90\silix_80\t" + row)
            if n != 0:
                
                o.write(f"{sample_name}\t{group_1}\t{group_2}\t{group_2}\t{group_2}\t{group_2}\t{group_2}\t" + row)
