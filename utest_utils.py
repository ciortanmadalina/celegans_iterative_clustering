import numpy as np
import pandas as pd
from scipy import sparse, io
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from collections import Counter
from IPython.display import clear_output, Image, display
from collections import OrderedDict
from scipy import stats
import scipy.stats as stats

def utest_scores(df, clusters, genes, top = 5, alpha = 0.05):
    """
    This method returns arrays of u test values of the distribution of the most distinctive genes
    from the input genes array (pass df.columns for complete analysis)
    If Top != None returns only top 5 most distinctive genes
    """
    #Array keeping the cluster names (keys)
    clusterNames= list(zip(*sorted(Counter(clusters).items())))[0]

    rankings_gene_names = []
    rankings_gene_pvals = []
    # iterate through clusters and calculate utest current cluster vs rest
    for igroup in clusterNames:
        idx = np.where(clusters == igroup)[0]
        idx_rest = np.where(clusters != igroup)[0]
        pvals = []
        cols = genes
        for gene in genes:
            if gene in df.columns:
                data = np.array(df[gene].tolist())
                stat, p = stats.mannwhitneyu(data[idx], data[idx_rest])
                pvals.append(p)
            
        # Return only first top scores, not all combinations
        if top!= None:
            # take only top genes most different, pval < alpha
            scores_sort_idx = np.argsort(pvals)[:top]
            pvals = np.array(pvals)[scores_sort_idx]
            cols = cols[scores_sort_idx]
            
        rankings_gene_names.append(cols)
        rankings_gene_pvals.append(pvals)
    return np.array(rankings_gene_names), np.array(rankings_gene_pvals)



def plot_utest(clusters, rankings_gene_names, rankings_gene_pvals, figname = 'u.pdf'):
    #Array keeping the cluster names (keys)
    clusterNames= list(zip(*sorted(Counter(clusters).items())))[0]
    
    n_panels_x = 4
    n_panels_y = np.ceil(len(clusterNames) / n_panels_x).astype(int)
    # n_panels_y = np.ceil(20 / n_panels_x).astype(int)
    print(n_panels_x, n_panels_y)

    fig = plt.figure(figsize = (n_panels_x * 4, n_panels_y * 4))
    left = 0.2/n_panels_x
    bottom = 0.13/n_panels_y
    gs = gridspec.GridSpec(nrows=n_panels_y,
                           ncols=n_panels_x,
                           left=left,
                           right=1-(n_panels_x-1)*left-0.01/n_panels_x,
                           bottom=bottom,
                           top=1-(n_panels_y-1)*bottom-0.1/n_panels_y,
                           wspace=0.22,
                           hspace=0.4)

    ax0 = None
    ymin = np.Inf
    ymax = -np.Inf
    for count, group_name in enumerate(clusterNames):
        ax = fig.add_subplot(gs[count])
        gene_names = rankings_gene_names[count]
        scores = rankings_gene_pvals[count]

        plt.grid()
        plt.scatter(np.arange(len(gene_names)),scores, s=30, alpha = 0.5)
        plt.xticks(np.arange(len(gene_names)), gene_names, rotation = 20)


        ax.set_title('Cluster {} vs. rest'.format(group_name))
        if count >= n_panels_x * (n_panels_y - 1):
            ax.set_xlabel('top distinctive genes')

        # print the 'score' label only on the first panel per row.
        if count % n_panels_x == 0:
            ax.set_ylabel('p val u score')
        ymin = np.min(scores)
        ymax = np.max(scores)
        ymax += 0.3*(np.max(scores)-np.min(scores))
        ax.set_ylim(ymin, ymax)
        ns = list(zip(*sorted(Counter(clusters).items())))[1]
        ns = np.array(ns)
        fig.suptitle(
            f'Total clusters {len(ns)}, valid(#>=10) {len(ns[ns>=10])}, invalid {len(ns[ns<10])}, outliers {np.sum(ns[ns<10])}',
            size=12 )
        plt.savefig(f'report/{figname}')
    
def utest_ficher_score(df, clusters, genes):
    """
    This method returns the fischer p value for the statistical test comparing 
    the expression of given genes in each cluster vs the rest.
    """
    #Array keeping the cluster names (keys)
    clusterNames= list(zip(*sorted(Counter(clusters).items())))[0]

    ficher_pvals = []
    # iterate through clusters and calculate utest current cluster vs rest
    for igroup in clusterNames:
        idx = np.where(clusters == igroup)[0]
        idx_rest = np.where(clusters != igroup)[0]
        pvals = []
        for gene in genes:
            if gene in df.columns:
                data = np.array(df[gene].tolist())
                stat, p = stats.mannwhitneyu(data[idx], data[idx_rest], alternative='greater')
                pvals.append(p)
        _, fischer_p = stats.combine_pvalues(pvals)

        ficher_pvals.append(fischer_p)
    return np.array(ficher_pvals)
        
def to_csv(clusters, df_matrix, name):
    result = df_matrix.copy()
    result ['cluster'] = clusters
    result = result[['cluster']]
    result.to_csv(f'output/{name}')