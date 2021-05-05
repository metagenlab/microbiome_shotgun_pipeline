import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import patches as mpatches
import matplotlib.colors as mcolor
from matplotlib.colors import LogNorm


def map_group_to_colors(group, palette):
    """
    Create a dictionary that maps a grouping variable to a color list
    """
    sorted_list = sorted(list(set(group)))
    rgb = sns.color_palette(palette, len(sorted_list)).copy()
    hex_colors = [mcolor.to_hex(color) for color in rgb]
    return dict(zip(sorted_list, hex_colors))


rank = snakemake.params.rank
method = snakemake.params.clustering
dist = snakemake.params.distance
tb = pd.read_csv(snakemake.input[0], sep='\t')

# get read counts matrix with samples as lines
subset = tb[tb[f'{rank}_taxid'] != 'na']
matrix = pd.pivot_table(subset, index=['sample', 'bodysite', 'type'], columns=f'{rank}',
                        values='read_counts', fill_value=0, aggfunc='sum').reset_index()
matrix = matrix.sort_values(by=list(matrix.columns), ascending=False)
metadata = matrix[['sample', 'bodysite', 'type']]
metadata.to_csv(snakemake.output.metadata, sep=',', index=None)
matrix.drop(['sample', 'bodysite', 'type'], axis=1).to_csv(snakemake.output.read_counts, sep=',', index=None)

# get read counts matrix with samples as column, get grouping variables colors
m = pd.pivot_table(subset, index=f'{rank}', columns='sample', values='read_counts', fill_value=0, aggfunc='sum')
m = m.sort_values(by=list(m.columns), ascending=False)

# get color table
samples2bodysite = dict(zip(metadata['sample'], metadata['bodysite']))
samples2type = dict(zip(metadata['sample'], metadata['type']))
body2colors = map_group_to_colors(metadata['bodysite'], 'hls')
type2colors = map_group_to_colors(metadata['type'], 'Paired')
colors = {}
for n, sample in enumerate(samples2bodysite):
    colors[n] = {}
    colors[n]['sample'] = sample
    colors[n]['type'] = type2colors[samples2type[sample]]
    colors[n]['bodysite'] = body2colors[samples2bodysite[sample]]
colors = pd.DataFrame.from_dict(colors, orient='index')
colors = colors.set_index('sample')

# plot clustermap
patches_bodysite = [mpatches.Patch(color=color, label=label) for label, color in body2colors.items()]
patches_type = [mpatches.Patch(color=color, label=label) for label, color in type2colors.items()]
g = sns.clustermap(m, norm=LogNorm(), figsize=(11.7, 8.27), row_cluster=False, col_colors=colors,
                   method=f'{method}', metric=f'{dist}', cmap='OrRd')
g.ax_heatmap.legend(handles=patches_bodysite, loc=2, title='bodysite', bbox_to_anchor=(-0.3, 0.8))
g.ax_col_colors.legend(handles=patches_type, loc=2, title='type', bbox_to_anchor=(-0.3, -0.2))
plt.title(loc='left', label=f'hierarchically-clustered heatmap, method={method}, metric={dist}')
g.savefig(snakemake.output.plot, dpi=400)
