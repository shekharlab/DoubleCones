from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
import samap.utils
from samalg import SAM
import pandas as pd
import scanpy as sc
import numpy as np
import seaborn as sn 
import matplotlib.pyplot as plt
from scipy.stats import zscore
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from sklearn.manifold import MDS 
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from functools import reduce
import os
import glob
import anndata as ad
from multiprocessing import Pool
from itertools import repeat

def MeanTable(results):
    # element-wise mean
    plt.rcParams['figure.figsize'] = [12, 10]
    tables = [result[1] for result in results]
    mean = pd.DataFrame(np.mean(tables, axis=0), columns=tables[0].columns, index=tables[0].index)
    PlotHeatmap(mean)
    
def StdTable(results):
    # element-wise std
    plt.rcParams['figure.figsize'] = [12, 10]
    tables = [result[1] for result in results]
    mean = pd.DataFrame(np.std(tables, axis=0), columns=tables[0].columns, index=tables[0].index)
    PlotHeatmap(mean)

# the score
def ScoreAlignment(df):

    df_nan = df.copy(deep = True)

    # Make interspecies comparisons NA
    for row in df_nan:
        for column in df_nan:
            if row[:2] == column[:2]:
                df_nan.loc[row,column] = np.nan

    # Rearraging Data
#     new_order = [x for x in new_order if x in df_nan.columns.values]

    # Reorder the rows and columns
#     df_reordered = df_nan.reindex(new_order)[new_order]
    df_reordered = df_nan

    # Metrics calculation
    goodsum_ = []
    badsum_ = []
    for row in df_reordered.index:
        for column in df_reordered.columns:
            if row[3:] == column[3:]:
                value1 = df_reordered.loc[row, column]
                if not pd.isna(value1):  # Check if value is not NaN
                    goodsum_.append(value1)  # Add value to the sum
            else:
                value2 = df_reordered.loc[row, column]
                if not pd.isna(value2):  
                    badsum_.append(value2) 

    min_ = -1*len(badsum_)
    max_ = 1*len(goodsum_)
    # score = sum(goodsum_) - sum(badsum_)
    
    # New metric ranges from -1 to 1 where 0 means similar enrichment of both good and bad edge weights
    accuracy = np.mean(goodsum_) - np.mean(badsum_) # old score: accuracy = (score - min_)/(max_ - min_)
    
    return accuracy, sum(goodsum_)/max_, -sum(badsum_)/min_

def unique(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

def SAMapTrial(filenames, gnnm = None, plot = False, remove_genes = None, keys = None, **kwargs):

    # Load objects
    objects = [ad.read_h5ad(file) for species, file in filenames.items()]
    objects
    
    # Remove genes
    if remove_genes is not None:
        objects2 = []
        for adata in objects:
            other_genes = [name for name in adata.var_names if not name in remove_genes]
            objects2.append(adata[:, other_genes])
        objects = objects2.copy()
    
    # Take a sample
    dsList = [DownsampleAnndata(object, **kwargs) for object in objects]
    
    # Run SAM
    samList = [SamPreprocessing(object) for object in dsList]
    samapDict = dict((identlist[i], samList[i]) for i, species in enumerate(specieslist))
    samapDict
    
    # Build SAMap object
    sm = RunSAM(samapDict, gnnm = gnnm)
    
    # Run SAMap
    sm, MappingTable, gpf = RunSAMap(sm, 
                                     NUMITERS = 3, 
                                     neigh_from_keys = dict((species, True) for species in identlist))

    if plot:
        plt.rcParams['figure.figsize'] = [3, 3]
        PlotUMAP(sm)
        
        plt.rcParams['figure.figsize'] = [12, 10]
        PlotHeatmap(MappingTable, sm)
        
    return sm, MappingTable, gpf

def SamPreprocessing(adata):
    sam=SAM(counts = adata)
    sam.preprocess_data(
                        sum_norm="cell_median",
                        norm="log",
                        thresh_low=0.0,
                        thresh_high=0.96, # Key parameter (this was set to 0.96 by default SAMap)
                        min_expression=1,
                    )
    sam.run(
                        preprocessing="StandardScaler",
                        npcs=100,
                        weight_PCs=False,
                        k=20,
                        n_genes=3000,
                        weight_mode='rms',
                        projection="none" # turn off umap calculation
                    )
    return sam

def DownsampleAnndata(adata, downsample = 100, cluster_key = "annotated", seed = 1):
    adatas = [adata[adata.obs[cluster_key] == clust] for clust in adata.obs[cluster_key].astype('category').cat.categories]

    for dat in adatas:
        if dat.n_obs > downsample:
             sc.pp.subsample(dat, n_obs=downsample, random_state=seed)

    adata_downsampled = adatas[0].concatenate(*adatas[1:])
    
    return adata_downsampled

def SamapBootstrap(filenames, gnnm, save = None):
    
    # delete preprocessed files if they exist
    for ident, file in filenames.items(): 
        file_to_delete = file.split('.h5ad')[0]+'_pr.h5ad'
        if os.path.isfile(file_to_delete):
            os.remove(file_to_delete)
            print('Removed ' + file_to_delete)
    
    sm = RunSAM(filenames.copy(), 
                gnnm = gnnm)
    
    sm, MappingTable, gpf = RunSAMap(sm, 
                                     NUMITERS = 3, 
                                     neigh_from_keys = dict((species, True) for species in filenames.keys()))#{'ze':True, 'ch':True, 'li':True, 'op':True, 'rn':True, 'hs':True})

    if save is not None:
        samap.utils.save_samap(sm, save)

    plt.rcParams['figure.figsize'] = [12, 10]
    PlotHeatmap(MappingTable, sm)
#     gene_pairs = gpf.find_all(align_thr=0.10)
#     print(gene_pairs)
    
    return sm, MappingTable, gpf

def RunSAM(filenames, keys = None, blast_dir = 'maps/', gnnm = None):

    if keys is None:
        sm = SAMAP(
                    filenames,
                    f_maps = blast_dir,
                    keys = dict((species, 'annotated') for species in filenames.keys()),
                    gnnm = gnnm, 
                    save_processed=False #if False, do not save the processed results to `*_pr.h5ad`
                )
    else: 
        sm = SAMAP(
                    filenames,
                    f_maps = blast_dir,
                    names = keys,
                    keys = dict((species, 'annotated') for species in filenames.keys()),
                    save_processed=False #if False, do not save the processed results to `*_pr.h5ad`
                )
    return sm

def RunSAM2(filenames, keys = None, blast_dir = 'maps/', gnnm = None):

    sams = []
    for name, file in filenames.items():
        sam=SAM()
        print("Loading ", file)
        sam.load_data(file)
        sam.preprocess_data(
                            sum_norm="cell_median",
                            norm="log",
                            thresh_low=0.0,
                            thresh_high=1, # Key parameter (this was set to 0.96 by default SAMap)
                            min_expression=1,
                        )
        sam.run(
                            preprocessing="StandardScaler",
                            npcs=100,
                            weight_PCs=False,
                            k=20,
                            n_genes=3000,
                            weight_mode='rms'
                        )
        sams.append(sam)

    print(dict(zip(filenames.keys(), sams)))
    
    sm = SAMAP(
            dict(zip(filenames.keys(), sams)),
            f_maps = blast_dir,
            gnnm = gnnm, 
#             names = keys,
            keys = dict((species, 'annotated') for species in filenames.keys())
        )

    return sm

def RunSAMap(sm, n_top = 0, **kwargs):
    
    # Run SAMap
    sm.run(**kwargs) #pairwise=pairwise, NUMITERS = NUMITERS)

    # Mapping
    keys = dict((species, 'annotated') for species in sm.ids)
    D,MappingTable = get_mapping_scores(sm, keys, n_top = n_top)
    
    gpf = GenePairFinder(sm,keys=keys)
    
    return sm, MappingTable, gpf

def MapTypes(filenames, keys = None, NUMITERS=3):
    
    # Run SAM
    sm = RunSAM(filenames, keys)
    
    # Run SAMap    
    return RunSAMap(sm, NUMITERS=NUMITERS)

def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64))
    return out

def RowNorm(dataframe):
    return dataframe.div(dataframe.sum(axis=1), axis=0)
    
def translate_feature_space(sm, ref, pred, average = True, scale = True, 
                            rownorm = False, scale_prediction = True, gnnm_refined = None):
    
    # Get gene names
    genes1 = [x for x in sm.gns if x.startswith(ref)]
    genes2 = [x for x in sm.gns if x.startswith(pred)]
    
    # Get average expression
    if average:
        X = grouped_obs_mean(sm.samap.adata[sm.samap.adata.obs['species'].isin([ref])][:, genes1], "annotated")
        Y = grouped_obs_mean(sm.samap.adata[sm.samap.adata.obs['species'].isin([pred])][:, genes2], "annotated")
    else: 
        X_adata = sm.samap.adata[sm.samap.adata.obs['species'].isin([ref])][:, genes1]
        X = pd.DataFrame(X_adata.X.toarray().T, index=genes1, columns=X_adata.obs.index)
        print(X.shape)
        Y_adata = sm.samap.adata[sm.samap.adata.obs['species'].isin([pred])][:, genes2]
        Y = pd.DataFrame(Y_adata.X.toarray().T, index=genes2, columns=Y_adata.obs.index)
        print(Y.shape)
    
    # Scale 
    if scale:
        X = pd.DataFrame(scaler.fit_transform(X.T).T, columns=X.columns, index=X.index)
#         X = X.apply(zscore, axis=1)
        X = X.fillna(0)
        Y = pd.DataFrame(scaler.fit_transform(Y.T).T, columns=Y.columns, index=Y.index)
#         Y = Y.apply(zscore, axis=1)
        Y = Y.fillna(0)
    
    # Generate prediction based on reference species and SAMap homology matrix
    if gnnm_refined is None:
        A = sm.gnnm_refined[np.ix_(np.flatnonzero(np.char.startswith(sm.gns, pred)),
                                   np.flatnonzero(np.char.startswith(sm.gns, ref)))]
    else:
        A = gnnm_refined[np.ix_(np.flatnonzero(np.char.startswith(sm.gns, pred)),
                                np.flatnonzero(np.char.startswith(sm.gns, ref)))]
        
    homology_matrix = pd.DataFrame(A.toarray(), 
                                   columns=genes1, 
                                   index=genes2)
    
    # Normalize homology matrix so that rows sum to 1 or cols sum to 1
    if(rownorm): 
        homology_matrix = RowNorm(homology_matrix)
    else: 
        homology_matrix = RowNorm(homology_matrix.T).T
    
    homology_matrix = homology_matrix.fillna(0)
    
    prediction = pd.DataFrame(homology_matrix @ X, 
                              index=genes2)
    
    # Scale prediction
    if(scale_prediction):
        prediction = pd.DataFrame(scaler.fit_transform(prediction.T).T, 
                                  columns=prediction.columns, 
                                  index=prediction.index)
        prediction = prediction.fillna(0)

    # Check that X is in right order
    if X.index.to_list() == genes1:
        print("The lists are identical")
    else:
        print("The lists are not identical")
    
    return X, Y, prediction, homology_matrix

def standardize_expression_matrices(pd_dict, ids, use_intersection = True):
    
    
    for k,v in enumerate(pd_dict):
        
        # Copy
        pd_dict[k] = pd_dict[k].copy(deep=True)
        
        # Add species name to each df
        pd_dict[k].columns = ids[k] + '_' + pd_dict[k].columns
        
        # Remove genes with zero variance
        pd_dict[k] = pd_dict[k].loc[~(pd_dict[k]==0).all(axis=1)]

    # Merge by intersection
    if use_intersection:
        merged = reduce(lambda x, y: pd.merge(x, y, how = 'inner', left_index=True, right_index=True), pd_dict)
    else:
        merged = pd.concat(pd_dict, axis=1)
        merged = merged.fillna(0)
    
    # Print
    print(merged.head())
    print('Shape: ', merged.shape)
    
    return merged

def linear_combination(sm, ref, pred, scale = False, normalize = True, use_A = False, test = False):
    
    # Get gene names
    genes1 = [x for x in sm.gns if x.startswith(ref)]
    genes2 = [x for x in sm.gns if x.startswith(pred)]
    
    # Get average expression
    X = grouped_obs_mean(sm.samap.adata[sm.samap.adata.obs['species'].isin([ref])][:, genes1], "annotated")
    Y = grouped_obs_mean(sm.samap.adata[sm.samap.adata.obs['species'].isin([pred])][:, genes2], "annotated")
    
    if scale:
        X = X.apply(zscore, axis=1)
        X = X.fillna(0)
        Y = Y.apply(zscore, axis=1)
        Y = Y.fillna(0)
    
    # Generate prediction based on reference species and SAMap homology matrix
    A = sm.gnnm_refined[np.ix_(np.flatnonzero(np.char.startswith(sm.gns, pred)),
                               np.flatnonzero(np.char.startswith(sm.gns, ref)))]
    homology_matrix = pd.DataFrame(A.toarray(), columns=genes1, index=genes2)
    
    # Check that X is in right order
    if X.index.to_list() == genes1:
        print("The lists are identical")
    else:
        print("The lists are not identical")
    
    if use_A:
        prediction = pd.DataFrame(A @ X, index=genes2, columns=X.columns)
    else: 
        prediction = pd.DataFrame(homology_matrix @ X, index=genes2)

    if test: 
        return prediction
    
    prediction['bias'] = 1
    
    # Apply linear regression
    w = np.linalg.inv(np.transpose(prediction) @ (prediction)) @ (np.transpose(prediction) @ Y)
    return prediction, Y
    #row_names = dict((idx, ele) for idx,ele in enumerate(w.columns))
    row_names = dict((idx, ele) for idx,ele in enumerate(X.columns.values))
    w = w.rename(index = row_names)
    
    # Normalization (each set of weights will sum to 1)
    if normalize:
#         w[w.columns] = w[w.columns] / w[w.columns].sum()
        wnorm=(w-w.min())/(w.max()-w.min()) # Min max normalization
    
    # Heatmap
    ax = sn.heatmap(data=wnorm, 
                   annot=True, 
                   annot_kws={"fontsize":8}, 
                   cmap='coolwarm',
                   center=0)

    plt.xlabel(pred) # x-axis label with fontsize 15
    plt.ylabel(ref) # y-axis label with fontsize 15
    plt.show()
    
    return X, Y, w, prediction, homology_matrix

def GoodnessOfFit(w, prediction, Y, pred = "Inferred", ref = "True", rsquared = True):
    best = prediction @ w.to_numpy() 

    # Compute R^2
    if rsquared: 
        R_squared = [[pearsonr(best.iloc[:,i], Y.iloc[:,j]).statistic**2 for i in range(len(best.columns))] for j in range(len(Y.columns))]
    else: 
        R_squared = [[mean_squared_error(Y.iloc[:,j], best.iloc[:,i], squared=False) for i in range(len(best.columns))] for j in range(len(Y.columns))]
    R_squared_matrix = pd.DataFrame(np.array(R_squared))
    
    # Annotate rows and columns
    R_squared_matrix.columns = w.columns
    R_squared_matrix.index = w.columns
    
    # Heatmap
    ax = sn.heatmap(data=R_squared_matrix, 
               annot=True, 
               annot_kws={"fontsize":8},
               cmap='coolwarm',
              center = np.min(R_squared_matrix))
    
    plt.xlabel(pred) # x-axis label with fontsize 15
    plt.ylabel(ref) # y-axis label with fontsize 15
    
    return R_squared_matrix

def GeneSimHeatmap(gpf, row_normalize = False):

    gene_pairs = gpf.find_all(align_thr=0)
    list_similarities = gene_pairs.loc[:, ~gene_pairs.columns.str.endswith(('_pval1', '_pval2'))].count()

    gene_summary = pd.DataFrame({'type1': [i.split(';', 1)[0] for i in list_similarities.index], 
                                  'type2': [i.split(';', 1)[1] for i in list_similarities.index],
                                  'n_genes':list_similarities})

    # make unique, sorted, common index
    idx = sorted(set(gene_summary['type1']).union(gene_summary['type2']))

    # reshape
    gene_matrix = (gene_summary.pivot(index='type1', columns='type2', values='n_genes')
       .reindex(index=idx, columns=idx)
       .fillna(0, downcast='infer')
       .pipe(lambda x: x+x.values.T)
     )

    # Row normalize
    if row_normalize:
        row_sums = gene_matrix.sum(axis=1)
        gene_matrix = gene_matrix / row_sums[:, np.newaxis]

    sn.heatmap(data=gene_matrix, 
               annot=True, 
               annot_kws={"fontsize":6}) 
    
    return gene_matrix


def PlotPCA(data, color = True, binary = None):
    # From https://builtin.com/machine-learning/pca-in-python
    # For subplots https://stackoverflow.com/questions/20073017/return-a-subplot-from-a-function

    X = data.copy(deep = True)
    if binary is not None:
        X[data >= binary] = 1
        X[data < binary] = 0
    
    # Standardizing the features
    x = StandardScaler().fit_transform(X.T)

    # PCA
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents
                 , columns = ['principal component 1', 'principal component 2'])
    principalDf.shape
    celltypes = X.columns.values.tolist()
    principalDf['celltype'] = celltypes
    finalDf = principalDf
    print(finalDf.head())
    # finalDf = pd.concat([principalDf, pd.DataFrame(target = X.columns.values.tolist())], axis = 1)
    
    if color:
        colors = [string.split('_')[1] for string in celltypes]
        for index, item in enumerate(colors):
            if item == "UV":
                colors[index] = "purple"
            if item == "rod":
                colors[index] = "grey"
            if item == "accessory":
                colors[index] = "cyan"
            if item == "principle":
                colors[index] = "magenta"
            if item == "principal":
                colors[index] = "magenta"
    else: 
        colors = "black"
            
    # Plot
    fig = plt.figure(figsize = (5,5))
    ax = fig.add_subplot(1,1,1) 
    ax.set_xlabel('PC1 (' + str(round(pca.explained_variance_ratio_[0]*100, 1)) + '%)', fontsize = 10)
    ax.set_ylabel('PC2 (' + str(round(pca.explained_variance_ratio_[1]*100, 1)) + '%)', fontsize = 10)
    ax.set_title('PCA', fontsize = 20)
#     ax.grid()
    ax.scatter(finalDf.loc[:, 'principal component 1'], 
               finalDf.loc[:, 'principal component 2'], 
               c = colors, 
               s = 50)
#     ax.legend(celltypes)

    for i, txt in enumerate(celltypes):
        ax.annotate(txt, (finalDf.loc[i, 'principal component 1'], finalDf.loc[i, 'principal component 2']))
        
    return finalDf
        

def PlotUMAP(sm, size=10):
    num_species = len(sm.ids)
    obs = sm.samap.adata.obs.replace('unassigned',np.NaN)
    annotation_column = "annotated;"*(num_species-1) + "annotated_mapping_scores"
    sm.samap.adata.obs['annotated'] = [i.split('_', 1)[1] for i in sm.samap.adata.obs[annotation_column].tolist()]
    sc.pl.umap(sm.samap.adata, 
           size=size,
           color=['species','annotated'],
           palette={'ze':'tab:green', 'op':'tab:blue','ch':'tab:red','li':'yellow','rn':'magenta', 'sq':'magenta','hs':'cyan',
                    'UV':'tab:purple','green':'tab:green','blue':'tab:blue','red':'tab:red','rod':'tab:grey', 
                    'OFF_bipolar':'black', 'ON_bipolar':'brown', 'AC':'black', 'RGC':'brown','HC':'brown', 
                    'SAC':'brown','VG3':'black',
                    'principle':'magenta', 'principal':'magenta', 'accessory':'cyan', 'opn1mw4/opn1lw1':'brown'})
    
def PlotHeatmap(MappingTable):
#     keys = dict((species, 'annotated') for species in sm.ids)
#     D,MappingTable = get_mapping_scores(sm, keys, n_top = 0)

    # Make interspecies comparisons NA
    for row in MappingTable:
        for column in MappingTable:
            if row[:2] == column[:2]:
                MappingTable.loc[row,column] = np.nan

    # Reorder
    ids = unique([column.split('_')[0] for column in MappingTable.columns])
    new_order = [species + '_' + type for type in ['UV', 'blue', 'green', 'opn1mw4/opn1lw1', 'red', 'principle','principal','accessory', 'rod'] for species in ids]
    new_order = [x for x in new_order if x in MappingTable.columns.values]
#     np.fill_diagonal(MappingTable.values, np.NaN)
    MappingTable = MappingTable.reindex(new_order)[new_order].round(2) # Round values
    sn.heatmap(data=MappingTable, 
               annot=True, 
               annot_kws={"fontsize":6})
    
def PlotMDS(MappingTable):
    mds = MDS(n_components=2, random_state=0) 
  
    # Fit the data to the MDS 
    # object and transform the data 
    X_transformed = pd.DataFrame(mds.fit_transform(MappingTable), columns = ['Comp1', 'Comp2'])
    X_transformed.index = MappingTable.index
    print(X_transformed)
    
    # Plot
    fig = plt.figure(figsize = (8,8))
    ax = fig.add_subplot(1,1,1)
    ax.grid()
    ax.set_xlabel('Component 1', fontsize = 15)
    ax.set_ylabel('Component 2', fontsize = 15)
    ax.set_title('MDS', fontsize = 20)
    
    colors = [string.split('_')[1] for string in MappingTable.index]
    for index, item in enumerate(colors):
        if item == "UV":
            colors[index] = "purple"
        if item == "rod":
            colors[index] = "grey"
        if item == "accessory":
            colors[index] = "cyan"
        if item == "principle":
            colors[index] = "magenta"
    
    ax.scatter(X_transformed.iloc[:,0], 
               X_transformed.iloc[:,1], 
               c = colors, 
               s = 50)
    
    for i, txt in enumerate(X_transformed.index):
        ax.annotate(txt, (X_transformed.loc[txt, 'Comp1'], X_transformed.loc[txt, 'Comp2']))

def PlotConnectivity(sm):
    c_matrix = pd.DataFrame(sm.samap.adata.obsp['connectivities'].A)
    cell_names = sm.samap.adata.obs.index

    # Annotations
    annotations = sm.samap.adata.obs['species'].astype(str) + '_' + sm.samap.adata.obs['annotated'].astype(str)
    new_order = [species + '_' + type for type in ['UV', 'blue', 'green', 'red', 'rod','principle','accessory','SAC','VG3','HC','AC'] for species in ['ze', 'ch', 'li', 'op', 'rn', 'hs']]
    new_order = [x for x in new_order if x in annotations.to_list()]
    annotations = pd.Categorical(annotations, categories=new_order)

    # Change the column names
    c_matrix.columns = annotations

    # Change the row indexes
    c_matrix.index = annotations

    # plt.rcParams['figure.figsize'] = [100, 100]
    plt.figure(figsize=(15, 15), dpi = 300) 

    # sn.heatmap(data=c_matrix.drop(columns=['annotation']), annot=False)
    ax = sn.heatmap(data=c_matrix.iloc[annotations.argsort(),annotations.argsort()], annot=False)
    line_indices = [annotations[annotations.argsort()].tolist().index(type) for type in new_order]
    ax.hlines(line_indices, *ax.get_xlim(), colors="white", linewidth=0.1)
    ax.vlines(line_indices, *ax.get_xlim(), colors="white", linewidth=0.1)
    ax