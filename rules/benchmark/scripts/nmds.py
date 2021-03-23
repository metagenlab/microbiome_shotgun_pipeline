import matplotlib.pyplot as plt
import pandas as pd
from sklearn import manifold
import seaborn as sns
import numpy as np
from scipy.spatial.distance import squareform
from scipy.spatial.distance import braycurtis
import itertools
from matplotlib.backends.backend_pdf import PdfPages

def Jaccard(X):
    """
    compute Jaccard similarity and then convert it into distance by subtracting it from 1.0.
    Args:
      X: input N x K data matrix. N ... the number of samples, K ... the number of features.
    Return:
      N x N data matrix. The value of (i,j) shows the distance between sample-i and sample-j.
    """
    X = np.array(X)
    n_samples = X.shape[0]
    n_distance = int(n_samples * (n_samples - 1) / 2)
    d_array = np.zeros(n_distance)
    for i, (idx1, idx2) in enumerate(itertools.combinations(range(n_samples), 2)):
        v1_nonzero_index = np.flatnonzero(X[idx1])
        v2_nonzero_index = np.flatnonzero(X[idx2])
        intersection = len(np.intersect1d(v1_nonzero_index, v2_nonzero_index))
        union = len(np.union1d(v1_nonzero_index, v2_nonzero_index))
        d_array[i] = 1.0 - (float(intersection) / float(union))
    return squareform(d_array)


def BrayCurtis(X):
    """
    compute Bray-Curtis dissimilarity.
    Args:
      X: input N x K data matrix. N ... the number of samples, K ... the number of features.
    Return:
      N x N data matrix. The value of (i,j) shows the distance between sample-i and sample-j.
    """
    X = np.array(X)
    n_samples = X.shape[0]
    n_distance = int(n_samples * (n_samples - 1) / 2)
    d_array = np.zeros(n_distance)
    for i, (idx1, idx2) in enumerate(itertools.combinations(range(n_samples),2)):
        d_array[i] = braycurtis(X[idx1], X[idx2])
    return squareform(d_array)

def plot_nmds(dist_mat, mat_index):
    """
    dist_mat: precomputed distance matrix
    mat_index: index names
    returns seaborn plot for the distance matrix
    """
    nmds = manifold.MDS(n_components=2, max_iter=3000, eps=1e-9, random_state=0,
                        dissimilarity="precomputed", n_jobs=1, metric=False)
    coords = nmds.fit(dist_mat).embedding_
    stress = round(nmds.stress_, 2)
    df = pd.DataFrame(coords, index=mat_index, columns=['NMDS1', 'NMDS2'])
    s = sns.scatterplot(data=df, x='NMDS1', y='NMDS2', hue='sample', style='sample', s=200)
    s.legend(loc=2, bbox_to_anchor=(1.05, 1), borderaxespad=0.1)
    l = min(df['NMDS1'])
    m = max(df['NMDS2'])
    s.text(l, m, f'2D stress: {stress}')

"""
plot heatmap
"""
all_samples = pd.read_csv(snakemake.input[0],sep='\t')
rank = snakemake.params.rank
threshold = snakemake.params.threshold
subset = all_samples[(all_samples[f'{rank}_taxid'] != 'na') & (all_samples['read_counts'] >= threshold)]
matrix = subset.pivot_table(index=f'{rank}', columns='sample', values='read_counts', fill_value=0, aggfunc='sum')
m = matrix.sort_values(by=list(matrix.columns), ascending=False)
plt.figure(figsize=(11.7, 8.27))
sns.set(font_scale=0.5)
plt.subplots_adjust(bottom=0.3, left=0.3, right=0.9)
hm = sns.heatmap(m, annot=True, fmt='d', cmap='OrRd')
hm.get_figure().savefig(snakemake.output.heatmap)
"""
plot NMDS
"""
x = subset.pivot_table(index='sample', columns=f'{rank}', values='read_counts', fill_value=0, aggfunc='sum')
tool = snakemake.wildcards.tool
with PdfPages(snakemake.output.nmds) as pdf:
    for distance in ['Jaccard', 'Bray-curtis']:
        plt.figure(figsize=(11.7, 8.27))
        plt.subplots_adjust(bottom=0.1, left=0.1, right=0.8)
        sns.set(font_scale=1.0)
        plt.title(f'tool: {tool}, distance: {distance}')
        if distance == 'Jaccard':
            distance_matrix = Jaccard(x)
        elif distance == 'Bray-curtis':
            distance_matrix = BrayCurtis(x)
        plot_nmds(distance_matrix, x.index)
        pdf.savefig()
        plt.close()
