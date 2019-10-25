
import pandas as pd
gold_standard=pd.read_csv(snakemake.input.gold_standard,sep='\t')
tool_output=pd.read_csv(snakemake.input.tool_out,sep='\t')
