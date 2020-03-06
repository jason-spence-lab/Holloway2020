import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
sc.settings.verbosity = 3			 # verbosity: errors (0), warnings (1), info (2), hints (3)
#sc.settings.set_figure_params(dpi_save=1200, dpi=1200)

greatestPalette = ["#f44336","#265bcf","#36b3ba","#ffeb3b","#e91e63","#00cc92","#4caf50","#ffb65c","#9c27b0","#03a9f4","#43d61f","#ff9800","#673ab7","#cddc39","#81bdc5","#ff5722"]
palestPalette = ["#fcc9c5","#acb4e2","#2effea","#fffbd6","#f7abc5","#b1dafb","#b5deb6","#ffe79e","#d88ae5","#90dbfe","#d5e9be","#ffd699","#bca6e3","#70eeff","#edf3ba","#ffccbd"]

#############################################################################
##         Flag to toggle rerunning the analysis or just plotting          ##
#############################################################################

rerun = True

#############################################################################
##                          Plot section flags                             ##
#############################################################################

do_ec_analysis = True
plot_raw = True
plot_normalized = False

redraw_featureplots = True
redraw_umaps = True

run_marker_analysis = True

expression_cutoff = 0.01  # Sets the threshold for min expression to be shown in plots

#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
mistorage_mount_point = 'Z:/'

#############################################################################
## Change this to list the sampleIDs of the samples you want to see        ##
#############################################################################
sample_list = ['2598-15','2598-10','2598-12','2598-2']


#############################################################################
## Change this to contain all genes you want to see expression for         ##
#############################################################################
#genes_of_interest = ['EPCAM', 'VIM', 'PDGFRA', 'ACTA2', 'TOP2A', 'STMN2', 'CDH5', 'KDR', 'CD34', 'TEK', 'NOTCH4', 'PLVAP', 'ESAM', 'GNG11', 'IGFBP4', 'MYL6', 'CALM1', 'RAMP2', 'ECSCR', 'TMSB10', 'FKBP1A', 'MEF2C', 'PDE5A', 'SFRP1', 'TSHZ2', 'EFNB2', 'GJA4', 'CSPG4', 'DES', 'RGS5', 'PDGFRB']

genes_of_interest = ['CDH5','KDR','EMCN','FLT1','TEK','CD34','ESAM','PLVAP','RAMP2','NOTCH4','VWF','PDPN','PROX1','ETV2','ERG',
					 'MEOX1','NKX2-3','FABP4','CA4','ADRB1','VIPR1','IGFBP5','IRX3','IRX5','CRHBP']

#############################################################################
##        Adjust these parameters to get a nicer looking UMAP plot         ##
#############################################################################
# UMAP arguments
num_neighbors_use = 50
num_pcs_use = 15
umap_spread = 1
umap_min_dist = 0.3
maxiter = 100
umap_gamma=1

# Louvain arguments
louv_res = 0.8

# PAGA arguments
size=20
paga_layout='fr'
threshold=0.005
node_size_scale=3

#############################################################################
##               DO NOT CHANGE ANYTHING BEYOND THIS POINT                  ##
#############################################################################




## Location to output the anndata h5ad files
raw_data_file = ''.join(['./data/Data_', '_'.join(sample_list), '.scanpy.raw.h5ad'])  # the file that will store the raw combined data
results_file = ''.join(['./data/Data_', '_'.join(sample_list), '.processed.h5ad'])  # the file that will store the analysis results
filtered_data_file = ''.join(['./data/Data_', '_'.join(sample_list), '.scanpy.filtered.h5ad'])  # the file that will store the raw combined data
endo_results_file = ''.join(['./data/Data_', '_'.join(sample_list), '.endothelial.processed.h5ad'])  # the file that will store the analysis results
neuro_results_file = ''.join(['./data/Data_', '_'.join(sample_list), '.neuronal.processed.h5ad'])  # the file that will store the analysis results

## Define function to generate a color gradient from a defined starting and ending color
def make_cmap(colors, position=None, bit=False):
	'''
	make_cmap takes a list of tuples which contain RGB values. The RGB
	values may either be in 8-bit [0 to 255] (in which bit must be set to
	True when called) or arithmetic [0 to 1] (default). make_cmap returns
	a cmap with equally spaced colors.
	Arrange your tuples so that the first color is the lowest value for the
	colorbar and the last is the highest.
	position contains values from 0 to 1 to dictate the location of each color.
	'''
	import matplotlib as mpl
	import numpy as np
	bit_rgb = np.linspace(0,1,256)
	if position == None:
		position = np.linspace(0,1,len(colors))
	else:
		if len(position) != len(colors):
			sys.exit("position length must be the same as colors")
		elif position[0] != 0 or position[-1] != 1:
			sys.exit("position must start with 0 and end with 1")
	if bit:
		for i in range(len(colors)):
			colors[i] = (bit_rgb[colors[i][0]],
						 bit_rgb[colors[i][1]],
						 bit_rgb[colors[i][2]])
	cdict = {'red':[], 'green':[], 'blue':[]}
	for pos, color in zip(position, colors):
		cdict['red'].append((pos, color[0], color[0]))
		cdict['green'].append((pos, color[1], color[1]))
		cdict['blue'].append((pos, color[2], color[2]))

	cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
	return cmap

annotation_dict = dict()

for line in open(''.join([mistorage_mount_point, '01_RNAseq_RAW_Data/single_cell_meta_data_table.tsv']), 'r'):
	#print(line)
	elem = str.split(line.rstrip())
	#print(elem)
	if elem:
		if elem[0] not in annotation_dict:
			annotation_dict[elem[0]] = elem[1:]

def Create_Scanpy_Anndata(mistorage_mount_point, sampleID):
	metadata_list = annotation_dict[sampleID][1:]
	newAdata = sc.read_10x_h5(''.join([mistorage_mount_point, annotation_dict[sampleID][0]]), genome='hg19')
	## Set gene names to be unique since there seem to be duplicate names from Cellranger
	newAdata.var_names_make_unique()
	## Add metadata for each sample to the observation (cells) annotations in the Anndata objects
	print('\nAdding Metadata for sample',sampleID,'\n')
	for field in metadata_list:
		field_list = str.split(field, ':')
		meta_name = field_list[0]
		meta_value = field_list[1]
		newAdata.obs[meta_name] = meta_value
	return(newAdata)

# function to get unique values 
def unique(list1): 
  
    # intilize a null list 
    unique_list = [] 
      
    # traverse for all elements 
    for x in list1: 
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 
    return(unique_list)

## Create my custom palette for FeaturePlots and define a matlplotlib colormap object
feature_colors = [(230,230,230), (35,35,142), (255,127,0)]
position=[0, expression_cutoff, 1]
my_feature_cmap = make_cmap(feature_colors, position=position, bit=True)
dot_colors = [(230,230,230), (153,0,0), (255,145,0)]
my_dot_cmap = make_cmap(dot_colors, position=position, bit=True)



## Read the raw Cellranger filtered data matrices into new Anndata objects
def runBasicAnalysis():
	
	first_adata = False
	
	if Path(raw_data_file).is_file():
		print(''.join(['Data_', '_'.join(sample_list), '.scanpy.raw.h5ad']), 'found, using this existing raw data file\n')
		adata = sc.read_h5ad(raw_data_file)
	else:
		print('\nNo existing h5ad raw data file found, reading in 10x h5 data for each sample\n')
		for sample in sample_list:
			if not first_adata:
				adata = Create_Scanpy_Anndata(mistorage_mount_point, sample)
				first_adata = True
			else:
				adata = adata.concatenate(Create_Scanpy_Anndata(mistorage_mount_point, sample))
		
		## Make cell names unique by adding _1, _2, _3 sequentially to each duplicated 10x barcode/name
		adata.obs_names_make_unique()
		## Write the raw combined dataset to disk so you can skip combining samples next time
		print('\nSaving raw combined sample data to', raw_data_file, '\n')
		adata.write(raw_data_file)
		open('./data/Prefiltered_gene_list.txt', 'w').write('\n'.join(adata.var_names.tolist()))
	
	## Basic filtering to get rid of useless cells and unexpressed genes
	
	sc.pp.filter_cells(adata, min_genes=1000)
	sc.pp.filter_genes(adata, min_cells=10)
	
	print('\nDoing initial filtering...\nKeeping', len(adata.obs_names), 'cells and', len(adata.var_names), 'genes.\n')
	
	mito_genes = adata.var_names.str.startswith('MT-')
	# Calculate the percent of genes derived from mito vs genome
	# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
	adata.obs['percent_mito'] = np.sum(
		adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
	# add the total counts per cell as observations-annotation to adata
	adata.obs['n_counts'] = adata.X.sum(axis=1).A1
	
	sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
				 jitter=0.4, multi_panel=True, save = '_preFiltering_plot.pdf', show = False)
	
	## Actually do the filtering.
	
	adata = adata[adata.obs['n_genes'] > 1000, :]   # Keep cells with more than 1000 genes
	adata = adata[adata.obs['n_genes'] < 9000, :]   # Keep cells with less than 5000 genes to remove most doublets
	adata = adata[adata.obs['n_counts'] < 60000, :] # Keep cells with less than 15000 UMIs to catch a few remaining doublets
	adata = adata[adata.obs['percent_mito'] < 0.5, :]   # Keep cells with less than 0.1 mito/genomic gene ratio
	sc.pp.filter_genes(adata, min_cells=50)	# Refilter genes to get rid of genes that are only in a tiny number of cells
	
	print('\nDoing final filtering...\nKeeping', len(adata.obs_names), 'cells and', len(adata.var_names), 'genes.\n')
	
	open('./data/Final_filtered_gene_list.txt', 'w').write('\n'.join(adata.var_names.tolist()))
	
	sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
				 jitter=0.4, multi_panel=True, save = '_postFiltering_plot.pdf', show = False)
	
	## Normalize the expression matrix to 10,000 reads per cell, so that counts become comparable among cells.
	# This corrects for differences in sequencing depth between cells and samples
	
	#sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.normalize_total(adata)
	
	## Log transform the data.
	
	sc.pp.log1p(adata)
	
	## Set the .raw attribute of AnnData object to the logarithmized raw gene expression for later use in differential testing and visualizations of gene expression.
	# We need to do this because the expression matrix will be rescaled and centered which flattens expression too much for some purposes
	
	adata.write(filtered_data_file)
	adata.raw = adata
	
	## Identify highly-variable genes based on dispersion relative to expression level.
	
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=6, min_disp=0.2)
	
	## Filter the genes to remove non-variable genes since they are uninformative
	
	adata = adata[:, adata.var['highly_variable']]
	
	## Regress out effects of total reads per cell and the percentage of mitochondrial genes expressed.
	
	sc.pp.regress_out(adata, ['n_counts'])
	
	## Scale each gene to unit variance. Clip values exceeding standard deviation 10 to remove extreme outliers
	
	sc.pp.scale(adata, max_value=10)
	
	## Run PCA to compute the default number of components
	
	sc.tl.pca(adata, svd_solver='arpack')
	
	## Rank genes according to contributions to PCs.
	
	sc.pl.pca_loadings(adata, show=False, components=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], save='_PCA-loadings.pdf')
	
	## Draw the PCA elbow plot to determine which PCs to use
	sc.pl.pca_variance_ratio(adata, log=True, n_pcs = 100, save = '_elbowPlot.pdf', show = False)
	
	## Compute nearest-neighbors
	
	sc.pp.neighbors(adata, n_neighbors=num_neighbors_use, n_pcs=num_pcs_use)
	
	## Calculate cell clusters via Louvain algorithm
	
	sc.tl.louvain(adata, resolution = louv_res)
	
	## Run UMAP Dim reduction
	
	sc.tl.umap(adata, min_dist=umap_min_dist, maxiter=maxiter, spread=umap_spread, gamma=umap_gamma)
	
	## Save the final clustered and dim-reduced data file as h5ad
	print('\nSaving processed Anndata data to', results_file, '\n')
	adata.write(results_file)
	return(adata)


if rerun:
	adata = runBasicAnalysis()
	adata.raw = sc.read_h5ad(filtered_data_file)
	print('Rerunning filtering and normalization...')
elif not rerun:
	print('Rerun toggled to "off"...\nMoving directly to subclustering and plotting\n') 
	adata = sc.read_h5ad(results_file)
	adata.raw = sc.read_h5ad(filtered_data_file)

orig_adata = sc.read_h5ad(raw_data_file)


if do_ec_analysis or plot_raw or plot_normalized:
	## And for endothelial cells
	endo_cells = adata[adata.obs['louvain']=='14'].obs.index.tolist()
	endo_adata = orig_adata[orig_adata.obs.index.isin(endo_cells)]
	print('\nSaving processed endothelial cell Anndata data to', endo_results_file, '\n')
	endo_adata.write(endo_results_file)
	
	adata_endo = sc.read_h5ad(endo_results_file)
	sc.pp.filter_genes(adata_endo, min_cells=5)	# Refilter genes to get rid of genes that are only in a tiny number of cells
	print('\nEndothelial cluster...\nContains', len(adata_endo.obs_names), 'cells, expressing', len(adata_endo.var_names), 'genes.\n')
	
	day_0_endo_cells = len(adata_endo[adata_endo.obs['age']=='0'].obs_names)
	day_3_endo_cells = len(adata_endo[adata_endo.obs['age']=='3'].obs_names)
	day_7_endo_cells = len(adata_endo[adata_endo.obs['age']=='7'].obs_names)
	day_14_endo_cells = len(adata_endo[adata_endo.obs['age']=='14'].obs_names)
	
	day_0_total_cells = len(adata[adata.obs['age']=='0'].obs_names)
	day_3_total_cells = len(adata[adata.obs['age']=='3'].obs_names)
	day_7_total_cells = len(adata[adata.obs['age']=='7'].obs_names)
	day_14_total_cells = len(adata[adata.obs['age']=='14'].obs_names)
	
	day_0_ratio = day_0_endo_cells / day_0_total_cells * 100
	day_3_ratio = day_3_endo_cells / day_3_total_cells * 100
	day_7_ratio = day_7_endo_cells / day_7_total_cells * 100
	day_14_ratio = day_14_endo_cells / day_14_total_cells * 100
	
	print('\nEndothelial Cell Counts:\nDay 0 =', day_0_endo_cells, '\nDay 3 =', day_3_endo_cells, '\nDay 7 =', day_7_endo_cells, '\nDay 14 =', day_14_endo_cells, '\n')
	print('\nTotal Cell Counts:\nDay 0 =', day_0_total_cells, '\nDay 3 =', day_3_total_cells, '\nDay 7 =', day_7_total_cells, '\nDay 14 =', day_14_total_cells, '\n')
	print('\nEndothelial Cell Percentage:\nDay 0 =', day_0_ratio, '\nDay 3 =', day_3_ratio, '\nDay 7 =', day_7_ratio, '\nDay 14 =', day_14_ratio, '\n')
	
	cellCountsOut = open('./figures/Endothelial_cell_count_data.txt', 'w')
	
	lineout = '\t'.join(['dataType', 'Day-0', 'Day-3', 'Day-7', 'Day-14', '\nraw_counts', str(day_0_endo_cells), str(day_3_endo_cells), str(day_7_endo_cells), str(day_14_endo_cells), '\ntotal_counts', str(day_0_total_cells), str(day_3_total_cells), str(day_7_total_cells), str(day_14_total_cells), '\npercent_EC', 	str(day_0_ratio), str(day_3_ratio), str(day_7_ratio), str(day_14_ratio)])
	
	cellCountsOut.write(lineout)
	
	cellCountsOut.close()






################################################################################
##                        RAW CELL COUNTS PER TIMEPOINT                       ##
################################################################################

if plot_raw:
	################################################################################
	## Pie chart, where the slices will be ordered and plotted counter-clockwise: ##
	################################################################################
	labels = 'Day 0', 'Day 3', 'Day 7', 'Day 14'
	sizeLabels = (str(day_0_endo_cells), str(day_3_endo_cells), str(day_7_endo_cells), str(day_14_endo_cells))
	sizes = [day_0_endo_cells, day_3_endo_cells, day_7_endo_cells, day_14_endo_cells]
	total = sum(sizes)
	
	explode = (0.05, 0.05, 0.05, 0.05)  # only "explode" day 3 since it's the peak of EC population
	
	piefig1, (pieax1, pietable) = plt.subplots(nrows=2, gridspec_kw=dict(height_ratios=[3,1]))

	#piefig1.suptitle('Endothelial Cells - Per Timepoint')
	
	pietable.axis("off")
	
	wedges, texts, autotexts = pieax1.pie(sizes, explode=explode, autopct=lambda p: '{:.0f}'.format(p * total / 100), shadow=False, startangle=90, colors=greatestPalette, pctdistance=1.1)
	
	pieax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
	
	pieax1.legend(wedges, labels, title="Timepoint", loc="lower right", bbox_to_anchor=(1, 0), fontsize='small')
	
	#pieax1.set_title('Endothelial cells - per timepoint')
	
	col_labels = ['Day 0', 'Day 3', 'Day 7', 'Day 14']
	table_sizes = [[day_0_endo_cells, day_3_endo_cells, day_7_endo_cells, day_14_endo_cells]]		   
	
	the_pie_table = pietable.table(cellText=table_sizes,
	          rowLabels=['No. ECs'],
	          colLabels=col_labels,
	          colColours=greatestPalette,
	          loc='center',
	          cellLoc='center')
	
	the_pie_table.set_fontsize(8)
	
	piefig1.savefig('./figures/EndothelialCells_by_stage_piechart.raw-counts.pdf')
	
	
	# Bar Graph for cells per stage within the endothelial cell cluster
	
	plt.rcdefaults()
	barfig, (barax, bartable) = plt.subplots(nrows=2, gridspec_kw=dict(height_ratios=[3,1]))
	
	#barfig.suptitle('Endothelial Cells - Per Timepoint')
	
	bartable.axis("off")
	
	bars = ('Day 0', 'Day 3', 'Day 7', 'Day 14')
	y_pos = np.arange(len(bars))
	
	barax.bar(y_pos, sizes, align='center', color=greatestPalette, width=0.5)
	barax.set_xticks(y_pos)
	barax.set_xticklabels(bars)
	barax.set_ylabel('Number of Cells')
	#barax.set_title('Endothelial cells - per timepoint')
	
	col_labels = ['Day 0', 'Day 3', 'Day 7', 'Day 14']
	table_sizes = [[day_0_endo_cells, day_3_endo_cells, day_7_endo_cells, day_14_endo_cells]]
	
	the_bar_table = bartable.table(cellText=table_sizes,
	          rowLabels=['No. ECs'],
	          colLabels=col_labels,
	          colColours=greatestPalette,
	          loc='center',
	          cellLoc='center')
	
	the_bar_table.set_fontsize(8)
	
	barfig.savefig('./figures/EndothelialCells_by_stage_bargraph.raw-counts.pdf')

















################################################################################
##            CELL COUNTS NORMALIZED TO TOTAL CELLS PER TIMEPOINT             ##
################################################################################

if plot_normalized:
	################################################################################
	## Pie chart, where the slices will be ordered and plotted counter-clockwise: ##
	################################################################################
	labels = 'Day 0', 'Day 3', 'Day 7', 'Day 14'
	sizeLabels = (str(day_0_endo_cells), str(day_3_endo_cells), str(day_7_endo_cells), str(day_14_endo_cells))
	sizes = [day_0_ratio, day_3_ratio, day_7_ratio, day_14_ratio]
	total = sum(sizes)
	
	explode = (0.05, 0.05, 0.05, 0.05)  # only "explode" day 3 since it's the peak of EC population
	
	piefig1, (pieax1, pietable) = plt.subplots(nrows=2, gridspec_kw=dict(height_ratios=[3,1]))
	
	#piefig1.suptitle('Percent Endothelial Cells - Per Timepoint')
	
	pietable.axis("off")
	
	wedges, texts, autotexts = pieax1.pie(sizes, explode=explode, autopct=lambda p: '{:.2f}'.format(p * total / 100), shadow=False, startangle=90, colors=greatestPalette, pctdistance=1.15)

	pieax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
	
	pieax1.legend(wedges, labels, title="Timepoint", loc="lower right", bbox_to_anchor=(1, 0), fontsize='small')
	
	#pieax1.set_title('Endothelial cells - per timepoint')
	
	col_labels = ['Day 0', 'Day 3', 'Day 7', 'Day 14']
	table_sizes = [[day_0_endo_cells, day_3_endo_cells, day_7_endo_cells, day_14_endo_cells],
				   [day_0_total_cells, day_3_total_cells, day_7_total_cells, day_14_total_cells],
				   ['{:.2f}'.format(day_0_ratio), '{:.2f}'.format(day_3_ratio), '{:.2f}'.format(day_7_ratio), '{:.2f}'.format(day_14_ratio)]]
				   
	
	the_pie_table = pietable.table(cellText=table_sizes,
	          rowLabels=['No. ECs', 'Total Cells', 'Pct ECs'],
	          colLabels=col_labels,
	          colColours=greatestPalette,
	          loc='center',
	          cellLoc='center')
	
	the_pie_table.set_fontsize(8)
	
	piefig1.savefig('./figures/EndothelialCells_by_stage_piechart.normalized-counts.pdf')
	
	#######################################################################
	## Bar Graph for cells per stage within the endothelial cell cluster ##
	#######################################################################
	plt.rcdefaults()
	barfig, (barax, bartable) = plt.subplots(nrows=2, gridspec_kw=dict(height_ratios=[3,1]))
	
	#barfig.suptitle('Percent Endothelial Cells - Per Timepoint')
	
	bartable.axis("off")
	
	cell_percentages = [day_0_ratio, day_3_ratio, day_7_ratio, day_14_ratio]
	
	bars = ('Day 0', 'Day 3', 'Day 7', 'Day 14')
	y_pos = np.arange(len(bars))
	
	barax.bar(y_pos, cell_percentages, align='center', color=greatestPalette, width=0.5)
	barax.set_xticks(y_pos)
	barax.set_xticklabels(bars)
	barax.set_ylabel('Percentage of Total Cells')
	#barax.set_title('Endothelial cells - per timepoint')
	
	col_labels = ['Day 0', 'Day 3', 'Day 7', 'Day 14']
	table_sizes = [[day_0_endo_cells, day_3_endo_cells, day_7_endo_cells, day_14_endo_cells],
				   [day_0_total_cells, day_3_total_cells, day_7_total_cells, day_14_total_cells],
				   ['{:.2f}'.format(day_0_ratio), '{:.2f}'.format(day_3_ratio), '{:.2f}'.format(day_7_ratio), '{:.2f}'.format(day_14_ratio)]]
	
	the_bar_table = bartable.table(cellText=table_sizes,
	          rowLabels=['No. ECs', 'Total Cells', 'Pct ECs'],
	          colLabels=col_labels,
	          colColours=greatestPalette,
	          loc='center',
	          cellLoc='center')
	
	the_bar_table.set_fontsize(8)
	
	barfig.savefig('./figures/EndothelialCells_by_stage_bargraph.normalized-counts.pdf')











'''

endo_data_points = adata[adata.obs['louvain']=='14'].obsm['X_umap']

other_data_points = adata.obsm['X_umap']


sns.set(style="white")


#print(adata[adata.obs['louvain']=='14'].raw[:, 'KDR'].X)
#print(adata.obs['louvain']=='14')

louv_list = adata.obs['louvain'].values.tolist()
endo_cat_list = []
for item in louv_list:
	if item =='13':
		endo_cat_list.append('endo')
	elif item !='13':
		endo_cat_list.append('other')

adata.obs['is_endo'] = endo_cat_list

kdr_expression = adata.raw[:, 'KDR'].X
cdh5_expresion = adata.raw[:, 'CDH5'].X
col1a1_expresion = adata.raw[:, 'COL1A1'].X


dataset = pd.DataFrame({'umap_1':adata.obsm['X_umap'][:,0],'umap_2':adata.obsm['X_umap'][:,1], 'KDR':adata.raw[:, 'KDR'].X, 'CDH5':adata.raw[:, 'CDH5'].X, 'COL1A1':adata.raw[:, 'COL1A1'].X, 'RAMP2':adata.raw[:, 'RAMP2'].X, 'MYL6':adata.raw[:, 'MYL6'].X, 'VIM':adata.raw[:, 'VIM'].X, 'Celltype':adata.obs['is_endo'].values, 'louvain':adata.obs['louvain'].values, 'age':adata.obs['age'].values})

kdr_dataset = dataset[dataset['KDR'] > 0 ]

df = pd.DataFrame({'MEOX1':adata.raw[:, 'MEOX1'].X, 'EBF3':adata.raw[:, 'EBF3'].X, 'CA4':adata.raw[:, 'CA4'].X, 'ADRB1':adata.raw[:, 'ADRB1'].X, 'Celltype':adata.obs['is_endo'].values})

g = sns.pairplot(df, hue="Celltype")

g.savefig('seaborn_scattermatrix_test.pdf')


s_plot = sns.relplot(x="umap_1", y="umap_2", hue="Celltype", size="CDH5", sizes=(4, 200), alpha=.5, palette="muted", height=6, data=dataset)
s_plot.savefig('seaborn_scatter_test.pdf')

R_plot = sns.relplot(x="MYL6", y="RAMP2", hue="Celltype", size="CDH5", sizes=(4, 200), alpha=.5, palette="muted", height=6, data=dataset)
R_plot.savefig('seaborn_relplot_test.pdf')

v_plot = sns.catplot(x="age", y="KDR", hue="Celltype", kind="violin", split=True, data=kdr_dataset);
v_plot.savefig('seaborn_violin_test.pdf')

b_plot = sns.catplot(x="age", y="KDR", kind="boxen", data=kdr_dataset.sort_values("age"));
b_plot.savefig('seaborn_boxen_test.pdf')


#divide_by_zero = 2983567938456 / 0

#endo_expression = adata[adata.obs['louvain']=='13'].raw[:, 'KDR'].X




'''









if redraw_umaps:
	print('\nRedrawing the umap plots...\n---------------------------\n')
	sc.pl.umap(adata, color='louvain', save = '_clusterIdentity.pdf', show = False, legend_loc = 'on data', edges = True, edges_color = 'lightgrey', edges_width = 0.01, size = 20, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
	sc.pl.umap(adata, color='louvain', save = '_clusterIdentity_noEdge.pdf', show = False, legend_loc = 'on data', edges = False, edges_color = 'lightgrey', edges_width = 0.01, size = 20, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
	sc.pl.umap(adata, color=['louvain', 'age'], save = '_clusterIdentity_age.pdf', show = False, legend_loc = 'right margin', edges = False, edges_color = 'lightgrey', edges_width = 0.01, size = 20, palette = greatestPalette, alpha = 0.95, legend_fontsize=6)
	sc.pl.umap(adata, color='age', save = '_age.pdf', show = False, legend_loc = 'right margin', edges = False, size = 20, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color='tissue', save = '_tissue.pdf', show = False, legend_loc = 'right margin', edges = False, size = 20, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color='sex', save = '_sex.pdf', show = False, legend_loc = 'right margin', edges = False, size = 20, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color='gel', save = '_hydrogel.pdf', show = False, legend_loc = 'right margin', edges = False, size = 20, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color='media', save = '_media.pdf', show = False, legend_loc = 'right margin', edges = False, size = 20, palette = greatestPalette, alpha = 0.95)
	sc.pl.umap(adata, color='sampleName', save = '_sample.pdf', show = False, legend_loc = 'right margin', edges = False, size = 20, palette = greatestPalette, alpha = 0.95)
	
	## Run PAGA to get nicer cell starting positions for UMAP
	sc.tl.paga(adata, groups='louvain')
	sc.pl.paga(adata, color='louvain', save='_paga-init-pattern.pdf', show=False, threshold=threshold, node_size_scale=node_size_scale, node_size_power=0.9, layout=paga_layout)
	
	## Setup multi-panel UMAP figure layouts
	
	ageFigumap = plt.figure(dpi=80, figsize=(18,12))
	#ax1 = ageFigumap.add_subplot(2,4,1)
	ax2 = ageFigumap.add_subplot(2,3,1)
	ax3 = ageFigumap.add_subplot(2,3,2)
	ax4 = ageFigumap.add_subplot(2,3,3)
	ax5 = ageFigumap.add_subplot(2,3,4)
	#ax6 = ageFigumap.add_subplot(2,4,6)
	ax0 = ageFigumap.add_subplot(2,3,5)
	
	sc.pl.umap(adata, color='louvain', show = False, legend_loc = 'on data', edges = False, size = 25, palette = greatestPalette, alpha = 0.95, legend_fontsize=8, ax=ax0)
	
	good_xlims = ax0.get_xlim()
	good_ylims = ax0.get_ylim() 
	
	#sc.pl.umap(adata[adata.obs['age']=='-5'], color='louvain', show = False, legend_loc = 'on data', edges = False, size = 20, alpha = 0.95, legend_fontsize=8, ax=ax1)
	sc.pl.umap(adata[adata.obs['age']=='0'], color='louvain', show = False, legend_loc = 'on data', edges = False, size = 20, alpha = 0.95, legend_fontsize=8, ax=ax2)
	sc.pl.umap(adata[adata.obs['age']=='3'], color='louvain', show = False, legend_loc = 'on data', edges = False, size = 20, alpha = 0.95, legend_fontsize=8, ax=ax3)
	sc.pl.umap(adata[adata.obs['age']=='7'], color='louvain', show = False, legend_loc = 'on data', edges = False, size = 20, alpha = 0.95, legend_fontsize=8, ax=ax4)
	sc.pl.umap(adata[adata.obs['age']=='14'], color='louvain', show = False, legend_loc = 'on data', edges = False, size = 20, alpha = 0.95, legend_fontsize=8, ax=ax5)
	#sc.pl.umap(adata[adata.obs['age']=='28'], color='louvain', show = False, legend_loc = 'on data', edges = False, size = 20, alpha = 0.95, legend_fontsize=8, ax=ax6)
	
	#ax1.set_xlim(good_xlims)
	#ax1.set_ylim(good_ylims)
	ax2.set_xlim(good_xlims)
	ax2.set_ylim(good_ylims)
	ax3.set_xlim(good_xlims)
	ax3.set_ylim(good_ylims)
	ax4.set_xlim(good_xlims)
	ax4.set_ylim(good_ylims)
	ax5.set_xlim(good_xlims)
	ax5.set_ylim(good_ylims)
	#ax6.set_xlim(good_xlims)
	#ax6.set_ylim(good_ylims)
	
	ax0.set_title('Louvain clusters: All ages')
	#ax1.set_title('Louvain clusters: Definitive Endoderm')
	ax2.set_title('Louvain clusters: Sphereoids')
	ax3.set_title('Louvain clusters: Day 3 HIOs')
	ax4.set_title('Louvain clusters: Day 7 HIOs')
	ax5.set_title('Louvain clusters: Day 14 HIOs')
	#ax6.set_title('Louvain clusters: Day 28 HIOs')
	
	ageFigumap.savefig('./figures/UMAP_clusterIdentity_age_panels.pdf')
	
	## Save the final clustered and dim-reduced data file as h5ad
	print('\nSaving processed Anndata with UMAP embeddings to', results_file, '\n')
	adata.write(results_file)







if redraw_featureplots:
	print('\nRedrawing the feature plots...\n---------------------------\n')
	## Do FeaturePlots for select genes	
	expressed_dict = dict()
	
	for gene in adata.raw.var_names.values.tolist():
		if gene not in expressed_dict:
			expressed_dict[str(gene)] = 1
			#print(gene)
	
	genes_to_plot = []
	
	for gene in genes_of_interest:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_featureplots.png', show = False, cmap = my_feature_cmap, size = 25, use_raw = True)
	
	genes_to_plot = []
	standard_marker_genes = ['EPCAM', 'VIM', 'ACTA2', 'TOP2A', 'WT1', 'CDH5', 'DCN', 'KDR', 'ITGAM', 'PTPRC', 'HBB', 'STMN2', 'S100B']
	
	for gene in standard_marker_genes:
		if gene in expressed_dict:
			genes_to_plot.append(gene)
		else:
			print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')
	
	print('Plotting standard marker genes:', ' '.join(genes_to_plot),'\n')
	
	sc.pl.umap(adata, color=genes_to_plot, save = '_standard_markers_featureplots.png', show = False, cmap = my_feature_cmap, size = 25, use_raw = True)


print('Starting gene expression plotting...\n------------------------------------\n')


print('Checking for expression of genes of interest\n')
expressed_dict = dict()

for gene in adata.raw.var_names.values.tolist():
	if gene not in expressed_dict:
		expressed_dict[str(gene)] = 1

genes_to_plot = []
	
for gene in genes_of_interest:
	if gene in expressed_dict:
		genes_to_plot.append(gene)
	else:
		print('Sorry,', gene, 'Is not expressed in this dataset or is invariable.\n')

print('Found cells expressing', ' '.join(genes_to_plot), '\n')

#sc.pl.violin(adata, keys=genes_to_plot, groupby='louvain', stripplot=False, save='_generalMarkers_plot.pdf', show = False, use_raw=True)
			 
#sc.pl.violin(adata, keys=['KRT8', 'COL1A1', 'CDH5', 'PLVAP'], groupby='louvain', stripplot=False, save='_generalMarkers-v2_plot.pdf', show = False, use_raw=True)

print('Setting up figures for violin and robo plots\n')

strip_rows = len(genes_to_plot) * 2.25

violinplotFig = plt.figure(dpi=120, figsize=(6,strip_rows))

axs = []

df = pd.DataFrame({'umap_1':adata.obsm['X_umap'][:,0],'umap_2':adata.obsm['X_umap'][:,1], 'louvain':adata.obs['louvain'].values})

for i in range(len(genes_to_plot)):
	ax = violinplotFig.add_subplot(len(genes_to_plot),1,i+1)
	axs.append(ax)
	df[genes_to_plot[i]] = adata.raw[:, genes_to_plot[i]].X
	sns.violinplot(x="louvain", y=genes_to_plot[i], data=df, ax=ax, palette=greatestPalette, scale='width', bw='silverman', linewidth=1, inner='quartile')
	ax.set_xlabel(' ')
	#sns.catplot(x="louvain", y=genes_to_plot[i], kind="boxen", data=df.sort_values("louvain"), ax=ax, palette=greatestPalette)

print('Saving violin plot\n')


#violinplotFig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
violinplotFig.tight_layout()
violinplotFig.savefig('./figures/Violinplot_GOI_by_cluster.pdf')

roboplotFig = plt.figure(dpi=120, figsize=(6,strip_rows))

axs = []

for i in range(len(genes_to_plot)):
	ax = roboplotFig.add_subplot(len(genes_to_plot),1,i+1)
	axs.append(ax)
	df[genes_to_plot[i]] = adata.raw[:, genes_to_plot[i]].X
	#sns.stripplot(x="louvain", y=genes_to_plot[i], data=df, jitter=True, linewidth=1, ax=ax, palette=greatestPalette)
	sns.catplot(x="louvain", y=genes_to_plot[i], kind="boxen", data=df.sort_values("louvain"), ax=ax, palette=greatestPalette)
	ax.set_xlabel(' ')

print('Saving robo plot\n')

roboplotFig.tight_layout()
roboplotFig.savefig('./figures/Roboplot_GOI_by_cluster.pdf')

boxplotFig = plt.figure(dpi=120, figsize=(6,strip_rows))

flierprops = dict(marker='2', markerfacecolor='none', markersize=4, linestyle='none', alpha=0.4)

for i in range(len(genes_to_plot)):
	ax = boxplotFig.add_subplot(len(genes_to_plot),1,i+1)
	axs.append(ax)
	df[genes_to_plot[i]] = adata.raw[:, genes_to_plot[i]].X
	#sns.stripplot(x="louvain", y=genes_to_plot[i], data=df, jitter=True, linewidth=1, ax=ax, palette=greatestPalette)
	sns.catplot(x="louvain", y=genes_to_plot[i], kind="box", data=df.sort_values("louvain"), ax=ax, palette=greatestPalette, linewidth=1, notch=False, flierprops=flierprops)
	ax.set_xlabel(' ')

print('Saving box plot\n')

boxplotFig.tight_layout()
boxplotFig.savefig('./figures/Boxplot_GOI_by_cluster.pdf')

sc.pl.dotplot(adata, var_names=genes_to_plot, groupby='louvain', use_raw=True, log=True, expression_cutoff=expression_cutoff, color_map=my_dot_cmap, save = '_GOI_dotPlot.pdf', show = False)

if run_marker_analysis:
	print("\nAll done with general workflow... now finding marker genes.\n")
	## Find marker genes via Wilxocon test based on Louvain cluster assignment
	# Create a simple plot to show the top 25 most significant markers for each cluster
	
	sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
	
	sc.tl.filter_rank_genes_groups(adata, groupby='louvain', use_raw=True, log=True, key_added='rank_genes_groups_filtered', min_in_group_fraction=0.25, min_fold_change=2, max_out_group_fraction=0.5)
	
	sc.pl.rank_genes_groups(adata, key='rank_genes_groups_filtered', n_genes=10, sharey=False, save = '_markerPlots.pdf', show = False)
	sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups_filtered',  n_genes=10, save = '_markerDotPlots.pdf', show = False)

print('\nDone with entire script execution')










