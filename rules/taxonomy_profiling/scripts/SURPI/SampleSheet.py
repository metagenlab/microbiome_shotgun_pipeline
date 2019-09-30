import pandas as pd
tb=pd.read_csv(snakemake.input[0],sep='\t')
sample_names=list(tb["SampleName"])
barcode=list(range(1,len(sample_names)+1))
data={"Barcode":barcode,"sample_name":sample_names}
samplesheet=pd.DataFrame(data)
samplesheet.to_csv(snakemake.output[0],sep=',',index=None)