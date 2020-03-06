'''
BASIC SINGLE CELL ANALYSIS SCRIPT
by Josh Wu
4 June, 2019

Relies heavily on the Scanpy Python module developed by the Theis Lab
Read more about Scanpy at https://scanpy.readthedocs.io/en/latest/index.html

Contains analysis of HIO samples for Emily Holloway

In progress ---
Moving to encapsulate parameters and relevant functions using class sca_set()
'''
import sys
sys.path.insert(0,'C:/Users/Josh/Desktop/sca/tools')

from sca_run import *
import csv
#from tools.pipelines import *

figdir = './figures_030520/'
an_run = sca_run()
#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
an_run.storage_mount_point = 'Z:/'

## IDs of samples as represented in the metadata table
an_run.sample_list = ['3088-1','3088-2']

# ## List of interesting genes
an_run.add_gene_list(markers = ['EPCAM','COL1A1','COL1A2','PDGFRA','DCN','TCF21','LUM',
								'PITX1','CDH1','CDH5','KDR','DCN','EMCN','FLT1','VEGFA',
								'TEK','CD34','ESAM','PLVAP','RAMP2','NOTCH4','VWF','SOX9',
								'SP7','RUNX2','COL2A1','ACAN','SOX17','CDX2','CDH1','TNC',
								'DCN','ALPL','SHH','IHH','COL9A1','ETV2','NKX2-3',
								'ANO1','KIT','HAND2'],
					 label='basic_list')

an_run.add_gene_list(markers = ['CDH5','ESAM','ERG','KDR','CD34','EMCN','VWF','PROX1','PDGFRB',
								'LYVE1','SOX17','CDH1','EPCAM','TCF21','VIM','DCN','PDX1',
								'PDGFRA','CSPG4','RGS5','VIL1','CDX2','MEOX1','NKX2-3','SATB2',
								'FABP4','CA4','ADRB1','VIPR1','IGFBP5','CRHBP','IRX3','IRX5'],
					 label='emily_genes')

an_run.add_gene_list(markers = ['F3','DLL1','NPY','GPX3','NRG1','EGF'],
					 label = 'addtional_genes')

an_run.add_gene_list(markers = ['CDH5','CDH1','DCN','NKX2-3','MEOX1','FABP4','VIPR1','ADRB1','CRHBP','IGFBP5','IRX3','IRX5'],
					 label = 'young_int_list')

an_run.add_gene_list(markers = ['CDH1','VIL1','CDX2','PDX1','SATB2','VIM','DCN',
								'TCF21','PDGFRA','CDH5','KDR','ESAM','ERG','CD34','EMCN'],
					 label='new_list', groupby_positions = ['2','3','0','1','4','5'])#['1','2','0','3','4','5'])

an_run.add_gene_list(markers = ['CA4','VIPR1','ADRB1','MEOX1','FABP4','NKX2-3','CRHBP','IRX5','IRX3'],
					 label='cluster_8_list')

# an_run.add_gene_list(markers = ['FCGR2C','PCAT14','SYN2'],
# 					 label='score_gene_test_list')

## Parameters used to filter the data - Mainly used to get rid of bad cells
an_run.set_filter_params(min_cells = 0, # Filter out cells 
						 min_genes = 500, # Filter out cells with fewer genes to remove dead cells
						 max_genes = 8000, # Filter out cells with more genes to remove most doublets
						 max_counts = 50000, # Filter out cells with more UMIs to catch a few remaining doublets
						 max_mito = 0.15) # Filter out cells with high mitochondrial gene content

## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 15, # Size of the local neighborhood used for manifold approximation
						   n_pcs = 11, # Number of principle components to use in construction of neighborhood graph
						   spread = 1, # In combination with min_dist determines how clumped embedded points are
						   min_dist = 0.4, # Minimum distance between points on the umap graph
						   resolution = 0.25, # High resolution attempts to increases # of clusters identified
						   cell_score_lists = ['lung_osec','intestine_osec','kidney_osec'])

an_run.set_plot_params(size = 10,
					   umap_obs = ['louvain','sampleName'],
					   exp_grouping = ['louvain'],
					   umap_feature_color = 'blue_orange',
					   umap_categorical_color = ['#33EB33','#005200','#1A6DCC','#000766','#009F00','#ED8D93'])#'#EF8E95','#7F7E7E',

## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
an_run.pipe_basic(figdir)#,load_save='adata_save.p',new_save=False)

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code
# New analysis parameters for the subset of parameters
analysis_params_ext = dict(n_neighbors = 9,
						n_pcs = 10,
						spread = 1,
						min_dist = 0.4,
						resolution = 0.4,
						cell_score_lists = ['lung_osec','intestine_osec','kidney_osec'])

an_run.umap_cluster_color = 'default'
an_run.size=160
#an_run.gene_dict['new_list']['groupby_positions'] = None
umap_categorical_color='default'

# an_run.pipe_ext(analysis_params_ext, figdir=figdir, extracted=['5'], load_save='adata_save.p', preprocess=False)

# #an_run.pipe_ext(analysis_params_ext, figdir=figdir, extracted=['8'], label='cluster_8/', load_save='adata_save.p')

# run_save = pickle.load(open(''.join([figdir,'extracted/adata_save.p']),"rb"))
# adata = run_save.adata.copy()
# adata.obs.loc[:,['lung_osec','intestine_osec','kidney_osec']].to_csv(''.join([figdir,'extracted_cell_scores.csv']))

# figdir = './figures_112019/young_int/'
# an_run.sample_list = ['2598-31','2757-2']

# an_run.pipe_basic(figdir)#,load_save='adata_save.p')

# figdir = './figures_112019/FACS/'
# an_run.sample_list = ['3078-1']

# ## Parameters used to filter the data - Mainly used to get rid of bad cells
# an_run.set_filter_params(min_cells = 0, # Filter out cells 
# 						 min_genes = 750, # Filter out cells with fewer genes to remove dead cells
# 						 max_genes = 2000, # Filter out cells with more genes to remove most doublets
# 						 max_counts = 5000, # Filter out cells with more UMIs to catch a few remaining doublets
# 						 max_mito = 0.10) # Filter out cells with high mitochondrial gene content
# an_run.min_dist = 0.1
# an_run.min_dist = 0.5

# an_run.pipe_basic(figdir)#,load_save='adata_save.p')