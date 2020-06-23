###Extract table from text file
f=open(snakemake.input[0],'r')
lines=[]
for l in f:
    lines.append(l)
tab=open(snakemake.output[0],'w')
tab.writelines(lines[4:])
tab.close()