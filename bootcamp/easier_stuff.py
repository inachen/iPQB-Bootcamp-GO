# coding: utf-8

# This file is for those of you who learned Python over the summer (you did that, right?).
# In this file, I've put all of the nitty-gritty details of what makes this website work.

# Except it doesn't work, because you need to write all the functions!

# Some of these functions will just make the website easier to use. Some of them are
# important for the enrichment and clustering tasks that your teammates are working on.

# If you need any help, ask your team or a TA.


# (don't delete this but don't worry about it either)
import os # a built-in module, for dealing with filenames
from . import app # this is part of the website guts



# These are all the files you have to work with. Go open them in a text  editor so you can
# get a feel for what they look like, because you need to parse each one to turn on a
# piece of the website.

# A list of yeast genes, with standard names and short descriptions.
GENE_INFO = os.path.join(app.root_path, 'data', 'gene_info.txt')

# A file that maps from GOID to name, aspect (process/function/component), etc
GO_INFO = os.path.join(app.root_path, 'data', 'go_info.txt')

# A two-column file that maps GOID to yeast genes
GO_MEMBERSHIP = os.path.join(app.root_path, 'data', 'go_membership.txt')

# A many-columned file that contains experimental data (yeast microarrays). Each column
# (after the first) is a different experiment, and each row is a gene. The values are log2
# ratios versus untreated control.
EXPERIMENT_FILE = os.path.join(app.root_path, 'data', 'experiment_data.txt')


# return a list or dictionary that maps from the id of an experiment (an int: 0, 1, ..)
# to a list of (systematic name, fold-change value) tuples
# e.g. [[('YAL001C', -0.06), ('YAL002W', -0.3), ('YAL003W', -0.07), ... ],
#       [('YAL001C', -0.58), ('YAL002W', 0.23), ('YAL003W', -0.25), ... ],
#        ... ]
def experiment():
	with open(EXPERIMENT_FILE) as experimentfile:
		initial_size = 0
		next(experimentfile)
		for geneline in experimentfile:
			elements = geneline.split('\t')
			if initial_size == 0:
				expt_list = [list() for i in range(len(elements))]
				initial_size = len(elements)-1
			else:
				for i in range(0, initial_size):
					expt_list[i].append((elements[0], float(elements[i+1])))
		return expt_list
# 	f=open(EXPERIMENT_FILE,'rU')
# 	fLineList=f.readlines()
# 	f.close()
# 	i=0
# 	dict={}
# 	for lines in fLineList:
# 		if i==0:
# 			i+=1
# 			continue
# 		
# 		match=re.search(r'(Y\w+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+([\w\.-]+)\s+', lines)
# 		print 
# 		if match:
# 			dict[match.group(1)]=match.group(id)
# 	for k,v in dict.items():
# 		print k,v
# 	return 


# map from a gene's systematic name to its standard name
# e.g. gene_name('YGR188C') returns 'BUB1'
def gene_name(gene):
	with open(GENE_INFO) as gene_info:
		for geneline in gene_info:
			splitline = geneline.split()
			if splitline[0] == gene:
				return splitline[1].strip()


# map from a systematic name to some info about the gene (whatever you want),
# e.g  'YGR188C' -> 'Protein kinase involved in the cell cycle checkpoint into anaphase'
def gene_info(gene):
	with open(GENE_INFO) as gene_info:
		for geneline in gene_info:
			splitgene = geneline.split('\t')
			if splitgene[0] == gene:
				infogene = splitgene[2]
				return infogene


# map from a systematic name to a list of GOIDs that the gene is associated with
# e.g. 'YGR188C' -> ['GO:0005694', 'GO:0000775', 'GO:0000778', ... ]
def gene_to_go(gene):
	with open(GO_MEMBERSHIP) as GO_MEM:
		gofiles = []
		for go in GO_MEM:
			gofiles.append(go)
		GOIDs = []
		for i in range(len(gofiles)):
			test = gofiles[i].split()
			if test[0] == gene:
				GOIDs.append(test[1])
		return GOIDs


# map from one of the GO aspects (P, F, and C, for Process, Function, Component),
# to a list of all the GOIDs in that aspect
# e.g. 'C' -> ['GO:0005737', 'GO:0005761', 'GO:0005763', ... ]
def go_aspect(aspect):
	goids = []
	with open(GO_INFO) as go_info:
		for goline in go_info:
			splitgo = goline.split('\t')
			if splitgo[2] == aspect:
				goids.append(splitgo[0])
		return goids


# map from a GOID (e.g. GO:0005737) to a *tuple* of the term and term definition
# e.g. 'GO:0005737' -> ('cytoplasm', 'All of the contents of a cell... (etc)'
def go_info(goid):
	with open(GO_INFO) as go_info:
		for goline in go_info:
			splitgo = goline.split('\t')
			if splitgo[0] == goid:
				termtuple = (splitgo[1].strip(), splitgo[3].strip())
				break
		return termtuple


# the reverse of the gene_to_go function: map from a GOID
# to a list of genes (systematic names)
# e.g. 'GO:0005737' -> ['YAL001C', 'YAL002W', 'YAL003W', ... ]
def go_to_gene(goid):
	with open(GO_MEMBERSHIP) as GO_MEM:
		gofiles = []
		for go in GO_MEM:
			gofiles.append(go)
		genes = []
		for i in range(len(gofiles)):
			test = gofiles[i].split()
			if test[1] == goid:
				genes.append(test[0])
		return genes

experiment()
