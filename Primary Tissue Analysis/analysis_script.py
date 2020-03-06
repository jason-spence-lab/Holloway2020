##### SINGLE CELL ANALYSIS SCRIPT
# Josh Wu
# 3 May, 2019
# Contains analysis for Emily Holloway
# Analysis of Lung, Intestine and Kidney Data

import sys
sys.path.insert(0,'C:/Users/Josh/Desktop')
#from scanpy_spence import *
import scanpy_sca as sca
import scanpy as sc
from matplotlib import rcParams
import matplotlib.pyplot as plt
import copy
import _pickle as pickle

figdir = './figures_020720/'
sca_dict = dict()

[adata,sca_dict] = pickle.load(open(''.join([figdir,'adata_save.p']),"rb"))

#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
sca_dict.update(mistorage_mount_point = 'Z:/')


sca_dict.update(sample_list = ['2444-1','2288-4','2321-2', # Lung
							   '2757-1','2761-1', # Kidney
							   '2598-28', '2598-29', '2250-1', '2250-2', '2598-24', '2182-5','2292-1','2511-2']) #Intestine

sca_dict.update(gene_list = ['CDH5','KDR','FLT1','NOS3','VWF','EMCN','CDH1','KRT8','EPCAM','ITGAM',
							 'PTPRC','COL1A1','COL1A2','PDGFRA','S100B','STMN2','TUBB3','IGFBP5','ASIC2','ERG',
							 'PDPN','FLT4','THY1'])

sca_dict.update(organ_genes = ['CA4','ADRB1','VIPR1','MYOC','MEOX1','FABP4','EBF3','NKX2-3','CRHBP','ASIC2','IRX5','IRX3'])

sca_dict.update(dot_list = ['VEGFA','KDR','CDH5','EMCN','ESAM','NR2F2','EPHB4','EFNB2','DLL4','LYVE1','PROX1',
							'WNT2', 'FMO2', 'ABLIM3', 'INMT', 'FGFR4', 'RSPO2', 'NKX2-1', 'HHIP', 'CDH1', 'TCF21'])

sca_dict.update(gene_list_2 = ['CA4','ADRB1','VIPR1','MYOC','MEOX1','FABP4','EBF3','SOX17',
						'NKX2-3','CRHBP','ASIC2','IGFBP5','IRX5','IRX3','CDH5','GFAP'])

Kidney = dict(tissue = ['Kidney'],
			  min_genes = 500,
			  max_genes = 7000,
			  max_counts = 30000,
			  max_mito = 0.15)

Lung = dict(tissue =['Lung-Distal'],
			min_genes = 1000,
	    	max_genes = 5000,
			max_counts = 15000,
			max_mito = 0.10)

Intestine = dict(tissue =['Duodenum','Ileum'],
				min_genes = 1000,
	    		max_genes = 5000,
				max_counts = 15000,
				max_mito = 0.10)
organ_filters = dict(lung = Lung,
					 kid = Kidney,
					 int = Intestine)

sca_dict.update(param_dict=dict(organ_filters = organ_filters,
								min_cells = 0,
								n_neighbors = 15,
								n_pcs = 11,
								min_dist = 0.4,
								resolution = 0.25,
								cell_score_lists=['neuro_score','epithelium_score','mesenchymal_score','immune_score','endothelial_score']))#,'neuro_score']))

# [adata,sca_dict] = sca.load_data(sca_dict)
# adata = sca.filter_specific_genes(adata,text_file = 'problematic_blood_genes.txt')
# #sc.settings.set_figure_params(fontsize=12)
# #rcParams['figure.figsize'] = 4, 4

# [adata,sca_dict] = sca.filter_data(adata,sca_dict)
# [adata,sca_dict] = sca.preprocess_data(adata,sca_dict)
# #adata = adata[:,adata.var_names = ]
# pickle.dump([adata,sca_dict],open(''.join([figdir,'adata_save.p']),"wb"),protocol=4)

cell_score_keys = [''.join([score_name,'_raw_scaled']) for score_name in sca_dict['param_dict']['cell_score_lists']]
vmin = [adata.obs[scores].values.min() for scores in cell_score_keys]
print(vmin)
vmax = [adata.obs[scores].values.max() for scores in cell_score_keys]
print(vmax)
######### All Together ############
# adata = sca.run_analysis(adata,sca_dict)
# adata = sca.plot_sca(adata,sca_dict,figdir = figdir)

###### KIDNEY #######
sca_dict['param_dict'].update(n_neighbors = 15,
							n_pcs = 11,
							min_dist = 0.4,
							resolution = 0.25,
							vmin=vmin,vmax=vmax)

analysis_param_ext = dict(n_neighbors = 11,
						  n_pcs = 10,
						  min_dist=0.3,
						  resolution = 0.4,
						  vmin=vmin,vmax=vmax)

extracted_kid = ['5']
sca_dict.update(extracted = extracted_kid)
adata_kid = adata[(adata.obs['tissue'] == 'Kidney'),:].copy()
adata_kid = sca.pipeline_ext(adata_kid,sca_dict,analysis_param_ext,extracted_kid,figdir = ''.join([figdir,'figures_kidney']))

###### LUNG ######
sca_dict['param_dict'].update(n_neighbors = 15,
							n_pcs = 11,
							min_dist = 0.4,
							resolution = 0.5,
							vmin=vmin,vmax=vmax)

analysis_param_ext = dict(n_neighbors = 6,
						  n_pcs = 9,
						  min_dist=0.3,
						  resolution = 0.4,
						  vmin=vmin,vmax=vmax)

extracted_lung = ['7']
sca_dict.update(extracted = extracted_lung)
adata_lung = adata[(adata.obs['tissue'] == 'Lung-Distal'),:].copy()
adata_lung = sca.pipeline_ext(adata_lung,sca_dict,analysis_param_ext,extracted_lung,figdir = ''.join([figdir,'figures_lung']))

###### INTESTINE #######
sca_dict['param_dict'].update(n_neighbors = 14,
							n_pcs = 10,
							min_dist = 0.5,
							resolution = 0.2,
							vmin=vmin,vmax=vmax)

analysis_param_ext = dict(n_neighbors = 6,
						  n_pcs = 9,
						  min_dist=0.3,
						  resolution = 0.4,
						  vmin=vmin,vmax=vmax)
extracted_int = ['8']
sca_dict.update(extracted = extracted_int)
adata_int = adata[((adata.obs['tissue'] =='Duodenum') | (adata.obs['tissue'] == 'Ileum')),:].copy()
adata_int = sca.pipeline_ext(adata_int,sca_dict,analysis_param_ext,extracted_int,figdir = ''.join([figdir,'figures_intestine']))
print("End of File")


#extracted_kid = ['1','6','9']
#extracted_lung = ['3','9']
#extracted_int = ['6','7']
#adata = adata[adata.obs['louvain'].isin([extracted])]


###### SNS PLOT TESTING
# adata_kid.obs['Endothelial'] = [1 if cell else 0 for cell in adata_kid.obs['louvain'].isin(sca_dict['extracted'])]

# CA4_endo, CA4_other = [],[]
# for cell,is_endo in zip(adata_kid.raw[:,'CRHBP'].X, adata_kid.obs['louvain'].isin(sca_dict['extracted'])):
# 	(CA4_endo if is_endo else CA4_other).append(cell)
	
# # df = pd.DataFrame({'CA4_endo':CA4_endo,'CA4_other':CA4_other,'Endothelial':adata_kid.obs['Endothelial'].values,
# # 	'endo_label':['Endothelial' if label else 'Other' for label in list(map(bool,adata_kid.obs['Endothelial'].values))] })

# print(CA4_endo)
# print(CA4_other)
# #g = sns.lmplot(x='CA4',y='Endothelial',data=df,scatter_kws={"s":20},y_jitter=0.02)
# g = plt.figure()
# sns.kdeplot(CA4_endo,color='r')#,hist=False#,ax=axes[0,1])
# sns.kdeplot(CA4_other,color='b')
# g.savefig(''.join([figdir,'figures_kidney/','log_regress_test.png']))

	# group_positions = [(0,3),(4,5),(6,8),(9,11),(12,13),(14,15),(16,16),(17,17),(18,20),(21,21),(22,22),(23,25)]
	# group_labels = ['Podocyte','Cap Mesenchyme','Renal Interstitium','Collecting Duct','Proximal Tubule','Loop of Henle','Distal Convoluted Tubule',
	# 				'Glomerular Mesangium','Endothelial','Hematopoetic','Epithelial','Mesenchymal']