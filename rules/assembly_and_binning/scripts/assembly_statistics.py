import pandas 
import statistics
import re 

assembly_files = snakemake.input["depth"]
trimmomatic_output_files = snakemake.input["trimmed_reads_stats"]
all_samples = pandas.read_csv(snakemake.params[0], sep="\t", index_col=0)

def parse_trimmomatic(summary_file):
    with open(summary_file, "r") as f:
        # Input Read Pairs: 4399916 Both Surviving: 3566108 (81.05%) Forward Only Surviving: 634940 (14.43%) Reverse Only Surviving: 43611 (0.99%) Dropped: 155257 (3.53%)
        for row in f:
            if "Input Read Pairs" in row:
                # ['4399916', '3566108', '634940', '43611', '155257']
                s = re.findall(": (\d+)", row)
                n_pairs = s[0]
                surviving_pairs = s[1]
                surviving_forward_only = s[2]
                surviving_reverse_only = s[3]
                dropped = s[4]
    return n_pairs, surviving_pairs, surviving_forward_only, surviving_reverse_only, dropped


o = open(snakemake.output[0], 'w')

header = ["sample_name", 
          "group_1",
          "group_2", 
          "raw_read_pairs",
          "trimmed_read_pairs",
          "trimmed_read_singletons",
          "assembly_size",
          "mean_depth",
          "median_depth",
          "assembly_size_depth_0_5",
          "assembly_size_depth_5_10",
          "assembly_size_depth_10_20",
          "assembly_size_depth_20_30",
          "assembly_size_depth_30_40",
          "assembly_size_depth_40_50",
          "assembly_size_depth_50_100",
          "assembly_size_depth_100_more",
          "assembly_size_depth_0_5_percent",
          "assembly_size_depth_5_10_percent",
          "assembly_size_depth_10_20_percent",
          "assembly_size_depth_20_30_percent",
          "assembly_size_depth_30_40_percent",
          "assembly_size_depth_40_50_percent",
          "assembly_size_depth_50_100_percent",
          "assembly_size_depth_100_more_percent"]

o.write('\t'.join(header) + '\n')

sample2trimmo_data = {}
for trimmo in trimmomatic_output_files:
    sample_name = trimmo.split("/")[1]
    sample2trimmo_data[sample_name] = {}
    n_pairs, surviving_pairs, surviving_forward_only, surviving_reverse_only, dropped = parse_trimmomatic(trimmo)
    sample2trimmo_data[sample_name]["n_pairs"] = n_pairs
    sample2trimmo_data[sample_name]["surviving_pairs"] = surviving_pairs
    sample2trimmo_data[sample_name]["surviving_forward_only"] = surviving_forward_only
    sample2trimmo_data[sample_name]["surviving_reverse_only"] = surviving_reverse_only
    sample2trimmo_data[sample_name]["dropped"] = dropped
    sample2trimmo_data[sample_name]["singletons"] = int(surviving_forward_only) + int(surviving_reverse_only)
    


for assembly in assembly_files:
    print(assembly)
    sample_name = assembly.split("/")[1]

    group_1 = all_samples.loc[sample_name, "group_1"]
    group_2 = all_samples.loc[sample_name, "group_2"]

    t = pandas.read_csv(assembly, sep="\t", header=None)

    n_read_pairs = sample2trimmo_data[sample_name]["n_pairs"]
    n_surviving_pairs = sample2trimmo_data[sample_name]["surviving_pairs"]
    n_read_singletons = sample2trimmo_data[sample_name]["singletons"]

    assembly_size = len(t.iloc[:,0])
    assembly_size_depth_5 = len(t[t.iloc[:,2] <= 5])
    assembly_size_depth_10 = len(t[(t.iloc[:,2] > 5) & (t.iloc[:,2] <= 10)])
    assembly_size_depth_20 = len(t[(t.iloc[:,2] > 10) & (t.iloc[:,2] <= 20)])
    assembly_size_depth_30 = len(t[(t.iloc[:,2] > 20) & (t.iloc[:,2] <= 30)])
    assembly_size_depth_40 = len(t[(t.iloc[:,2] > 30) & (t.iloc[:,2] <= 40)])
    assembly_size_depth_50 = len(t[(t.iloc[:,2] > 40) & (t.iloc[:,2] <= 50)])
    assembly_size_depth_100 = len(t[(t.iloc[:,2] > 50) & (t.iloc[:,2] <= 100)])
    assembly_size_depth_100_more = len(t[t.iloc[:,2] > 100])
    assembly_size_depth_5_percent = ( assembly_size_depth_5 / float(assembly_size) ) * 100
    assembly_size_depth_10_percent = ( assembly_size_depth_10 / float(assembly_size) ) * 100
    assembly_size_depth_20_percent = ( assembly_size_depth_20 / float(assembly_size) ) * 100
    assembly_size_depth_30_percent = ( assembly_size_depth_30 / float(assembly_size) ) * 100
    assembly_size_depth_40_percent = ( assembly_size_depth_40 / float(assembly_size) ) * 100
    assembly_size_depth_50_percent = ( assembly_size_depth_50 / float(assembly_size) ) * 100
    assembly_size_depth_100_percent = ( assembly_size_depth_100 / float(assembly_size) ) * 100
    assembly_size_depth_100_more_percent = ( assembly_size_depth_100_more / float(assembly_size) ) * 100
    mean_depth = statistics.mean(t.iloc[:,2])
    median_depth = statistics.median(t.iloc[:,2])

    o.write(f"{sample_name}\t{group_1}\t{group_2}\t{n_read_pairs}\t{n_surviving_pairs}\t{n_read_singletons}\t{assembly_size}\t{mean_depth}\t{median_depth}\t"
            f"{assembly_size_depth_5}\t{assembly_size_depth_10}\t{assembly_size_depth_20}\t"
            f"{assembly_size_depth_30}\t{assembly_size_depth_40}\t{assembly_size_depth_50}\t{assembly_size_depth_100}\t{assembly_size_depth_100_more}\t"
            f"{assembly_size_depth_5_percent}\t{assembly_size_depth_10_percent}\t{assembly_size_depth_20_percent}\t"
            f"{assembly_size_depth_30_percent}\t{assembly_size_depth_40_percent}\t{assembly_size_depth_50_percent}\t"
            f"{assembly_size_depth_100_percent}\t{assembly_size_depth_100_more_percent}\n")
     
