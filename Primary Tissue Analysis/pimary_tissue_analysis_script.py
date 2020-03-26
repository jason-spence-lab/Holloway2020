##### FETAl-ADULT TRACHEA COMPARISON ANALYSIS SCRIPT
# Josh Wu
# 20 June, 2019
# Contains analysis of primary lung, kidney and intestine tissue for Emily Holloway
# Special analysis script

# import sys
# sys.path.insert(0,'C:/Users/Josh/Desktop/sca/tools')
from sca_tools.sca_run import *
# from sca_run import *
import copy

figdir = './figures_021820/'
an_run = sca_run()

#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
#sca_dict.update(storage_mount_point = 'Z:/')
an_run.storage_mount_point = 'Z:/'

## IDs of samples as represented in the metadata table
an_run.sample_list = ['2444-1','2288-4','2321-2', # Lung
					  '2757-1','2761-1', # Kidney
					  '2598-28', '2598-29', '2250-1', '2250-2', '2598-24', '2182-5','2292-1','2511-2'] #Intestine

# ## List of interesting genes
an_run.add_gene_list(markers=['CDH5','KDR','FLT1','NOS3','VWF','EMCN','CDH1','KRT8','EPCAM','ITGAM',
							 'PTPRC','COL1A1','COL1A2','PDGFRA','S100B','STMN2','TUBB3','IGFBP5','ASIC2','ERG',
							 'PDPN','FLT4','THY1'], label='gene_list')

an_run.add_gene_list(markers=['CA4','ADRB1','VIPR1','MYOC','MEOX1','FABP4','EBF3','NKX2-3','CRHBP','ASIC2','IRX5','IRX3'], 
					 label='organ_genes')

an_run.add_gene_list(markers=['VEGFA','KDR','CDH5','EMCN','ESAM','NR2F2','EPHB4','EFNB2','DLL4','LYVE1','PROX1',
							 'WNT2', 'FMO2', 'ABLIM3', 'INMT', 'FGFR4', 'RSPO2', 'NKX2-1', 'HHIP', 'CDH1', 'TCF21'],
					 label='dot_list')

an_run.add_gene_list(markers=['CA4','ADRB1','VIPR1','MYOC','MEOX1','FABP4','EBF3','SOX17',
							  'NKX2-3','CRHBP','ASIC2','IGFBP5','IRX5','IRX3','CDH5','GFAP'],label='gene_list_2')

an_run.add_gene_list(markers=['EPCAM','CDH1','KRT18','KRT8','FXYD3','PERP','CLDN6','POSTN','DCN',
							  'TCF21','COL1A2','COL3A1','COL1A1','TAGLN','ACTA2','FHL1','PTPRC',
							  'PTPRC','CD37','CORO1A','LCP1','CD53','LAPTM5','S100B','ELAVL4',
							  'TUBB2B','STMN2','ASCL1','NNAT','GRP','MPZ','CDH5','CLDN5','ESAM','KDR','FLT1'],
					 label='celltype_scoring_genes')

an_run.set_plot_params(size = 5,
					   umap_obs = ['louvain','sampleName'],
					   exp_grouping = ['louvain'],
					   umap_feature_color='blue_orange')

## Load all samples together
adata = an_run.load_data()
adata = an_run.filter_specific_genes(adata,text_file='problematic_blood_genes.txt')

## Split by tissue and filter separately
adata_kid = adata[adata.obs['tissue'].isin(['Kidney'])].copy()
adata_lung = adata[adata.obs['tissue'].isin(['Lung-Distal'])].copy()
adata_int = adata[adata.obs['tissue'].isin(['Duodenum','Ileum'])].copy()

## Kidney filtering
an_run.set_filter_params(min_cells = 0,
						 min_genes = 500, # Filter out cells with fewer genes to remove dead cells
						 max_genes = 7000, # Filter out cells with more genes to remove most doublets
						 max_counts = 30000, # Filter out cells with more UMIs to catch a few remaining doublets
						 max_mito = 0.15) # Filter out cells with high mitochondrial gene content

adata_kid = an_run.filter_data(adata_kid)
an_run_kid = copy.copy(an_run)
an_run_kid.adata_postFiltered = adata_kid.copy()

## Lung filtering
an_run.set_filter_params(min_cells = 0,
						 min_genes = 1000, # Filter out cells with fewer genes to remove dead cells
						 max_genes = 5000, # Filter out cells with more genes to remove most doublets
						 max_counts = 15000, # Filter out cells with more UMIs to catch a few remaining doublets
						 max_mito = 0.10) # Filter out cells with high mitochondrial gene content

adata_lung = an_run.filter_data(adata_lung)
an_run_lung = copy.copy(an_run)
an_run_lung.adata_postFiltered = adata_lung.copy()

##Intestine filtering
an_run.set_filter_params(min_cells = 0,
						 min_genes = 1000, # Filter out cells with fewer genes to remove dead cells
						 max_genes = 5000, # Filter out cells with more genes to remove most doublets
						 max_counts = 15000, # Filter out cells with more UMIs to catch a few remaining doublets
						 max_mito = 0.1) # Filter out cells with high mitochondrial gene content

adata_int = an_run.filter_data(adata_int)
an_run_int = copy.copy(an_run)
an_run_int.adata_postFiltered = adata_int.copy()


## Combine for preprocessing steps (log normalizing, scaling, etc.)
adata = adata_kid.concatenate(adata_lung,adata_int).copy()
an_run.cell_score_lists = ['neuro_score','epithelium_score','mesenchymal_score',
						   'immune_score','endothelial_score'] ## Marker lists for cell type scoring
adata = an_run.preprocess_data(adata)

cell_score_keys = an_run.cell_score_lists
vmin = [adata.obs[scores].values.min() for scores in cell_score_keys]
vmax = [adata.obs[scores].values.max() for scores in cell_score_keys]

an_run_kid.vmin_list=an_run_lung.vmin_list=an_run_int.vmin_list=vmin
an_run_kid.vmax_list=an_run_lung.vmax_list=an_run_int.vmax_list=vmax

## Separate again for further downstream analyses (UMAP clustering, etc.)
an_run_kid.adata = adata[(adata.obs['tissue'] == 'Kidney'),:].copy()
an_run_lung.adata = adata[(adata.obs['tissue'] == 'Lung-Distal'),:].copy()
an_run_int.adata = adata[((adata.obs['tissue'] =='Duodenum') | (adata.obs['tissue'] == 'Ileum')),:].copy()

os.makedirs(os.path.dirname(''.join([figdir,'figures_kidney/'])), exist_ok=True)
os.makedirs(os.path.dirname(''.join([figdir,'figures_lung/'])), exist_ok=True)
os.makedirs(os.path.dirname(''.join([figdir,'figures_intestine/'])), exist_ok=True)

pickle.dump(an_run_kid,open(''.join([figdir,'figures_kidney/','adata_save.p']),"wb"),protocol=4)
pickle.dump(an_run_lung,open(''.join([figdir,'figures_lung/','adata_save.p']),"wb"),protocol=4)
pickle.dump(an_run_int,open(''.join([figdir,'figures_intestine/','adata_save.p']),"wb"),protocol=4)


########################
###### KIDNEY #######
########################
## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 15, # Size of the local neighborhood used for manifold approximation
						   n_pcs = 11, # Number of principle components to use in construction of neighborhood graph
						   spread = 1, # In combination with min_dist determines how clumped embedded points are
						   min_dist = 0.4, # Minimum distance between points on the umap graph
						   resolution = 0.5, # High resolution attempts to increases # of clusters identified9
						   cell_score_lists = ['neuro_score','epithelium_score','mesenchymal_score','immune_score','endothelial_score'])

## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
an_run.pipe_basic(figdir=''.join([figdir,'figures_kidney/']),load_save='adata_save.p')

analysis_params_ext = dict(n_neighbors = 11,
						n_pcs = 10,
						spread = 1,
						min_dist = 0.3,
						resolution = 0.4,
						cell_score_lists = ['neuro_score','epithelium_score','mesenchymal_score','immune_score','endothelial_score'])

an_run.size=20
an_run.pipe_ext(analysis_params_ext,figdir=''.join([figdir,'figures_kidney/']),extracted=['6'],load_save='adata_save.p',preprocess=False)


########################
###### LUNG #######
########################
## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 15, # Size of the local neighborhood used for manifold approximation
						   n_pcs = 11, # Number of principle components to use in construction of neighborhood graph
						   spread = 1, # In combination with min_dist determines how clumped embedded points are
						   min_dist = 0.4, # Minimum distance between points on the umap graph
						   resolution = 0.5, # High resolution attempts to increases # of clusters identified9
						   cell_score_lists = ['neuro_score','epithelium_score','mesenchymal_score','immune_score','endothelial_score'])
an_run.size=5
## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
# an_run.pipe_basic(figdir=''.join([figdir,'figures_lung/']),load_save='adata_save.p')

analysis_params_ext = dict(n_neighbors = 6,
						n_pcs = 9,
						spread = 1,
						min_dist = 0.3,
						resolution = 0.4,
						cell_score_lists = ['neuro_score','epithelium_score','mesenchymal_score','immune_score','endothelial_score'])

an_run.size=20
# an_run.pipe_ext(analysis_params_ext,figdir=''.join([figdir,'figures_lung/']),extracted=['6'],load_save='adata_save.p',preprocess=False)


########################
###### INTESTINE #######
########################
## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 14, # Size of the local neighborhood used for manifold approximation
						   n_pcs = 10, # Number of principle components to use in construction of neighborhood graph
						   spread = 1, # In combination with min_dist determines how clumped embedded points are
						   min_dist = 0.5, # Minimum distance between points on the umap graph
						   resolution = 0.2, # High resolution attempts to increases # of clusters identified9
						   cell_score_lists = ['neuro_score','epithelium_score','mesenchymal_score','immune_score','endothelial_score'])
an_run.size=5
## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
# an_run.pipe_basic(figdir=''.join([figdir,'figures_intestine/']),load_save='adata_save.p')

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code
# Will make more robust
# New analysis parameters for the subset of parameters
analysis_params_ext = dict(n_neighbors = 6,
						n_pcs = 9,
						spread = 1,
						min_dist = 0.3,
						resolution = 0.4,
						cell_score_lists = ['neuro_score','epithelium_score','mesenchymal_score','immune_score','endothelial_score'])

an_run.size=20
# an_run.pipe_ext(analysis_params_ext,figdir=''.join([figdir,'figures_intestine/']),extracted=['8'],load_save='adata_save.p',preprocess=False)
