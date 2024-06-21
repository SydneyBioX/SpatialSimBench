import anndata as ad
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.preprocessing import LabelEncoder
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import squidpy as sq
import pandas as pd
from scipy.stats import ks_2samp

#%%
# transition matrix

def get_spatial_network(num_sample=None, spatial=None, radius=None, coord_type="grid", n_rings=2, set_diag=False):
    spatial_adata = ad.AnnData(np.empty((num_sample, 1), dtype="float32"))
    spatial_adata.obsm["spatial"] = spatial
    # sq.gr.spatial_neighbors(spatial_adata, n_rings=n_rings, coord_type=coord_type, n_neighs=n_neighs, radius=radius,set_diag =set_diag)
    sq.gr.spatial_neighbors(spatial_adata, n_rings=n_rings, coord_type=coord_type, radius=radius, set_diag=set_diag,
                            delaunay=True)
    sn = spatial_adata.obsp["spatial_connectivities"]

    return sn


def get_onehot_ct(init_assign=None):
    label_encoder = LabelEncoder()
    integer_encoded = label_encoder.fit_transform(init_assign)
    onehot_encoder = OneHotEncoder(sparse=False)
    integer_encoded = integer_encoded.reshape(len(integer_encoded), 1)
    onehot_ct = onehot_encoder.fit_transform(integer_encoded)
    return onehot_ct.astype(np.float32)


# @numba.jit("float32[:, ::1](float32[:, ::1], float32[:, ::1])")
def get_nb_freq(nb_count=None, onehot_ct=None):
    #     nb_freq = onehot_ct.T @ nb_count
    nb_freq = np.dot(onehot_ct.T, nb_count)
    res = nb_freq / nb_freq.sum(axis=1).reshape(onehot_ct.shape[1], -1)
    return res


def get_trans(adata=None, ct=None):
    sn = get_spatial_network(num_sample=adata.obs.shape[0],
                             spatial=adata.obsm["spatial"], coord_type="generic")
    onehot_ct = get_onehot_ct(init_assign=ct)
    nb_count = np.array(sn * onehot_ct, dtype=np.float32)
    target_trans = get_nb_freq(nb_count=nb_count, onehot_ct=onehot_ct)
    return target_trans


# build adata

def build_adata(num_sample=None, spatial_loc=None, cell_type=None):
    adata = ad.AnnData(np.empty((num_sample, 1), dtype="float32"))
    adata.obsm["spatial"] = spatial_loc
    # if cell type is int, then we need to transform it
    adata.obs["celltype"] = cell_type
    adata.obs["celltype"] = adata.obs["celltype"].astype('category')
    sq.gr.spatial_neighbors(adata, coord_type="generic", set_diag=False, delaunay=True)
    # neighborhood enrichment matrix
    sq.gr.nhood_enrichment(adata, cluster_key="celltype")
    # centrality scores matrix
    sq.gr.centrality_scores(adata, cluster_key="celltype")
    return adata


# compare plot

def compare_plot(adata_real, adata_sim, real, sim, file_path):
    warnings.simplefilter(action='ignore', category=FutureWarning)

    # transition matrix
    transition_matrix_real = get_trans(adata=adata_real, ct=real)
    transition_matrix_sim = get_trans(adata=adata_sim, ct=sim)

    error = np.linalg.norm(transition_matrix_sim - transition_matrix_real)
    transition_matrix_real_ds = transition_matrix_real.flatten()
    transition_matrix_sim_ds = transition_matrix_sim.flatten()
    ks_stat_error, p_value = ks_2samp(transition_matrix_real_ds, transition_matrix_sim_ds)
    
    # neihborhood enrichment matrix
    target_enrich_real = adata_real.uns["celltype_nhood_enrichment"]["zscore"]
    target_enrich_scale_real = target_enrich_real/np.max(target_enrich_real)
    target_enrich_sim = adata_sim.uns["celltype_nhood_enrichment"]["zscore"]
    target_enrich_scale_sim = target_enrich_sim/np.max(target_enrich_sim)

    error_enrich = np.linalg.norm(target_enrich_sim - target_enrich_real)
    error_enrich_scale = np.linalg.norm(target_enrich_scale_sim - target_enrich_scale_real)
    
    target_enrich_real_ds = target_enrich_real.flatten()
    target_enrich_sim_ds = target_enrich_sim.flatten()
    ks_enrich, p_value = ks_2samp(target_enrich_real_ds, target_enrich_sim_ds)

    # centrality score matrix
    real_central_real = np.array(adata_real.uns["celltype_centrality_scores"])
    real_central_sim = np.array(adata_sim.uns["celltype_centrality_scores"])
    
    real_central_real_ds = real_central_real.flatten()
    real_central_sim_ds = real_central_sim.flatten()
    ks_central, p_value = ks_2samp(real_central_real_ds, real_central_sim_ds)

    error_central = np.linalg.norm(np.array(real_central_sim) - np.array(real_central_real))
    error_degree_centrality = np.linalg.norm(np.array(real_central_sim)[:, 0] - np.array(real_central_real)[:, 0])
    error_average_clustering = np.linalg.norm(np.array(real_central_sim)[:, 1] - np.array(real_central_real)[:, 1])
    error_closeness_centrality = np.linalg.norm(np.array(real_central_sim)[:, 2] - np.array(real_central_real)[:, 2])

    # List of values
    # values = [error, error_enrich_scale, error_central, error_degree_centrality, error_average_clustering, error_closeness_centrality]
    values = [ks_stat_error, ks_enrich, ks_central]
    
    
    
    # Corresponding labels
    # labels = ['TM', 'NWE', 'CSM', 'CSMa', 'CSMb', 'CSMc']
    labels = ['TM', 'NWE', 'CSM']
    
    # save csv file
    results_df = pd.DataFrame({'Metric': labels, 'Value': values})
    csv_file_path = file_path.replace('.png', '_metricsNew.csv') 
    results_df.to_csv(csv_file_path, index=False)

    # Bar plot using seaborn for better control over aesthetics
    # palette = ["#313795", "#4575b4", "#74add1", "#abd9e9", "#fdae61", "#f46d43"]
    palette = ["#313795", "#4575b4", "#74add1"]
    sns.barplot(x=labels, y=values, palette=palette)

    ax = plt.gca()
    ax.set_facecolor('white')

    # Get the current figure and set the background to white
    fig = plt.gcf()
    fig.set_facecolor('white')

    # Adjusting the plot appearance with matplotlib
    plt.ylim(0, max(values)+0.5)
    plt.title('spatial level')
    plt.ylabel('Value')
    plt.xlabel('Metric')  # Since your R code had a 'Metric' on x-axis

    # Custom theme adjustments
    plt.grid(False)
    for spine in plt.gca().spines.values():
        spine.set_visible(True)
        spine.set_edgecolor('black')
        spine.set_linewidth(0.2)

    plt.savefig(file_path, format='png', dpi=300)

    plt.close()

    return plt


def final_spatial(df, file_path):
    spatial_coords = np.array(df[['x', 'y']].values.tolist())
    num_sample_test = len(spatial_coords)
    real_list = np.array(df['real'].values.tolist())
    sim_list = np.array(df['sim'].values.tolist())

    adata_real = build_adata(num_sample_test, spatial_coords, real_list)
    adata_sim = build_adata(num_sample_test, spatial_coords, sim_list)

    final_plot = compare_plot(adata_real, adata_sim, real_list, sim_list,file_path)

    return final_plot


# final_spatial(pd.read_csv("scDesign3_cluster_loc.csv"), 'output\\test\\spatial.png')

