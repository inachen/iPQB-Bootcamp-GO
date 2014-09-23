import numpy as np
from StringIO import StringIO
# coding: utf-8

# this file is for those of you with a little more coding experience. You need to calculate
# go-enrichment scores using some gene expression values.

# Feel free to add even more features if you like--the backend code is simple to figure
# out. But don't forget to help your teammates, and to figure out your perturbation!




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

def calculate_enrichment(gene_data, go_to_genes, n=100):

    sorted_genes = sorted(gene_data, key=lambda tup: tup[1], reverse=True)
    top_genes=sorted_genes[:n]
    top_gene_names= list(zip(*top_genes)[0])

    score_list = dict((k,0) for k in test_goid.keys())
    # print score_list

    for g in top_gene_names:
        if g in score_list:
            score_list[g] 




    positive_enrichment_scores = []
    negative_enrichment_scores = []

    return positive_enrichment_scores,negative_enrichment_scores


# You can make your website fancier by creating a figure for each experiment
# (it's up to your what this displays). Install the mpld3 module with pip:
#    pip install mpld3
# and read their website to get started: http://mpld3.github.io/index.html

# uncomment this line to import the module after you've installed it
# import mpld3

# input:
# - a list of [(gene name, gene value) ... ]
#
# output:
# - a dictionary, created by the mpld3 module from a figure you've made (see below)
#
def plot_experiment(gene_data):
    mpld3_dict = None

    # When you've made your plot, convert it with the mpld3 library like so:
    # mpld3_dict = mpld3.fig_to_dict(fig)

    return mpld3_dict

calculate_enrichment(test_genes, test_goid, n=3)

