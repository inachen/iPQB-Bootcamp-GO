import numpy as np
from StringIO import StringIO
from scipy.stats import hypergeom
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

# coding: utf-8

# this file is for those of you with a little more coding experience. You need to calculate
# go-enrichment scores using some gene expression values.

# Feel free to add even more features if you like--the backend code is simple to figure
# out. But don't forget to help your teammates, and to figure out your perturbation!


# cold shock, dtt, h2o2, heat shock, hypo-osmotic shock, nitrogen depletion

# inputs:
# - a list of [(gene name, gene value) ... ]
# - a dictionary of {GOID: [list of genes associated with this term]}
# - an int, N, to use as a parameter
#
# outputs:
# - a list of GOIDs, with associated enrichment scores,
#   testing for positive enrichment (sorted from most significant to least)
# - another list of GOIDs and enrichment scores,
#   testing for negative enrichment (again sorted)
#

test_genes = [('gene1', 1), ('gene2', 2), ('gene3', -1), ('gene4', 0), ('gene5', 1)]
test_goid = {'goid1': ['gene1', 'gene2'], 'goid2': ['gene3', 'gene4', 'gene5']}


# x = num genes in smaller list with goid
# M = total length of dataset
# N = num genes in total dataset with goid
# n = passed number
def calculate_enrichment(gene_data, go_to_genes, n=100):

    # get top n gene list
    top_sorted_genes = sorted(gene_data, key=lambda tup: tup[1], reverse=True)
    top_genes=top_sorted_genes[:n]
    top_gene_names= list(zip(*top_genes)[0])

    # get bottom n gene list
    bot_sorted_genes = sorted(gene_data, key=lambda tup: tup[1])
    bot_genes=bot_sorted_genes[:n]
    bot_gene_names= list(zip(*bot_genes)[0])

    tot_genes = len(gene_data)

    # [+ top + goid (k), total genes (M), top genes (n), goid genes (N), score ]
    # create score dictionary
    top_score_list = dict((k,[0, tot_genes, n, len(go_to_genes[k]), 0]) for k in go_to_genes.keys())
    bot_score_list = dict((k,[0, tot_genes, n, len(go_to_genes[k]), 0]) for k in go_to_genes.keys())

    # calculate top hits
    for g in top_gene_names:
        for goid in top_score_list:
            if g in go_to_genes[goid]:
                top_score_list[goid][0] += 1

    # calculate bottom hits
    for g in bot_gene_names:
        for goid in bot_score_list:
            if g in go_to_genes[goid]:
                bot_score_list[goid][0] += 1

    positive_enrichment_scores = []
    negative_enrichment_scores = []

    # calculate scores
    for goid in top_score_list:
        top_score_list[goid][4] = hypergeom.logsf(top_score_list[goid][0]-1, top_score_list[goid][1], top_score_list[goid][2], top_score_list[goid][3])
        positive_enrichment_scores.append((goid, top_score_list[goid][4]))
        bot_score_list[goid][4] = hypergeom.logsf(bot_score_list[goid][0]-1, bot_score_list[goid][1], bot_score_list[goid][2], bot_score_list[goid][3])
        negative_enrichment_scores.append((goid, bot_score_list[goid][4]))
            
    # 

	positive_enrichment_scores_sorted = sorted(positive_enrichment_scores, key=lambda tup: tup[1])
	negative_enrichment_scores_sorted = sorted(negative_enrichment_scores, key=lambda tup: tup[1])

    return positive_enrichment_scores_sorted,negative_enrichment_scores_sorted


# You can make your website fancier by creating a figure for each experiment
# (it's up to your what this displays). Install the mpld3 module with pip:
#    pip install mpld3
# and read their website to get started: http://mpld3.github.io/index.html

# uncomment this line to import the module after you've installed it
import mpld3
import operator

# input:
# - a list of [(gene name, gene value) ... ]
#
# output:
# - a dictionary, created by the mpld3 module from a figure you've made (see below)
#
def plot_experiment(gene_data):
    gene_data.sort(key=operator.itemgetter(1))
    filtered_genes = [gene_tuple[0] for gene_tuple in gene_data]
    filtered_values = [gene_tuple[1] for gene_tuple in gene_data]
    filtered_values.reverse()
    filtered_genes.reverse()
    
    fig, ax = plt.subplots(1)
    scatter = ax.scatter(range(len(gene_data)), filtered_values, edgecolor='none', facecolor='#ff69b4')
    
    tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=filtered_genes)
    mpld3.plugins.connect(fig, tooltip)

    # When you've made your plot, convert it with the mpld3 library like so:
    mpld3_dict = mpld3.fig_to_dict(fig)

    return mpld3_dict


