import numpy as np
import scanpy as sc
import json
from cmdstanpy import CmdStanModel
import scipy


# Create the data object
def create_gmm_data(small_infer):
    data = {}
    data['C'] = len(small_infer)
    data['G'] = len(small_infer[0])
    data['prior_mu'] = [.9109, 1.001, 1.1167]
    data['infer'] = small_infer

    GMM_data = 'scHPF_GMM.data.json'

    with open(GMM_data, "w") as outfile:
        json.dump(data, outfile)
    return GMM_data


def run_gmm(stan_model, infer):
    GMM_data = create_gmm_data(infer)
    GMM_stan = stan_model
    GMM_model = CmdStanModel(stan_file=GMM_stan)
    GMM_model.name
    GMM_model.stan_file
    GMM_model.exe_file
    GMM_model.code()
    GMM_fit = GMM_model.variational(data=GMM_data, output_dir='.',
                                    show_console=False,
                                    require_converged=False,
                                    tol_rel_obj=.001,
                                    iter=20000)
    return GMM_fit


def get_t_tests(ppd, data):
    deletion_d, amplification_d, nothing_d = [], [], []
    for d_point in data:
        if d_point < .95:
            deletion_d.append(d_point)
        elif d_point > 1.05:
            amplification_d.append(d_point)
        else:
            nothing_d.append(d_point)
    deletion_p, amplification_p, nothing_p = [], [], []
    for d_point in ppd:
        if d_point < .95:
            deletion_p.append(d_point)
        elif d_point > 1.05:
            amplification_p.append(d_point)
        else:
            nothing_p.append(d_point)
    return [scipy.stats.ttest_ind(deletion_d, deletion_p),
            scipy.stats.ttest_ind(nothing_p, nothing_d),
            scipy.stats.ttest_ind(amplification_p,
                                  amplification_d)]


def create_hpf_data(adata_subsampled, X, mixture_means):
    data = {}
    data['C'] = len(X)
    data['G'] = len(X[0])
    data['a_prime'] = 1
    std_1 = np.std(np.sum(adata_subsampled.layers['counts'], axis=1))
    data['b_prime'] = float(np.mean(np.sum(adata_subsampled.layers['counts'],
                                           axis=1))/std_1)
    data['a'] = .3
    data['c_prime'] = 1
    std_0 = np.std(np.sum(adata_subsampled.layers['counts'], axis=0))
    data['d_prime'] = float(np.mean(np.sum(adata_subsampled.layers['counts'],
                                           axis=0))/std_0)
    data['c'] = .3
    data['eta_sigma'] = .5
    data['k'] = 1
    data['X'] = X
    data['mixture_mu'] = mixture_means
    mod_scHPF_data = 'mod_scHPF.data.json'
    with open(mod_scHPF_data, "w") as outfile:
        json.dump(data, outfile)
    return mod_scHPF_data


def run_hpf(stan_model, adata_subsampled, X, mixture_means):
    mod_scHPF_data = create_hpf_data(adata_subsampled, X, mixture_means)
    mod_scHPF_stan = stan_model
    mod_scHPF_model = CmdStanModel(stan_file=mod_scHPF_stan)
    mod_scHPF_model.name
    mod_scHPF_model.stan_file
    mod_scHPF_model.exe_file
    mod_scHPF_model.code()
    mod_scHPF_fit = mod_scHPF_model.variational(data=mod_scHPF_data,
                                                output_dir='.',
                                                show_console=False,
                                                require_converged=False,
                                                tol_rel_obj=.001, iter=10000)
    return mod_scHPF_fit


def subsample(adata, key, index, n):
    '''
    adata: AnnData object to subsample from
    key: key in adata.obs that contains metadata of interest
    index: value of key of interest, i.e. keep if adata.obs[key]==index
    n: number of cells to select from cells for which adata.obs[key]==index
    '''
    keep = [i for i in range(adata.shape[0]) if adata.obs[key][i] == index]
    keep_subsampled = np.random.choice(keep, size=n, replace=False)
    return adata[keep_subsampled, :]


def format_data(adata, infer, n_genes, n_cells_per):
    keep = []
    for i in range(adata.shape[1]):
        if adata.var_names[i] in infer.columns:
            keep.append(i)
    adata = adata[:, keep]
    sc.pp.highly_variable_genes(adata, n_top_genes=n_genes, subset=True)
    clonal_ad = []
    for i in np.unique(adata.obs["inferCNV_clones"]):
        clonal_ad.append(subsample(adata, "inferCNV_clones", i, n_cells_per))
    adata_subsampled = clonal_ad[0].concatenate(clonal_ad[1:])
    infer_subsampled = np.zeros((adata_subsampled.shape))
    for i in range(adata_subsampled.shape[0]):
        for j in range(adata_subsampled.shape[1]):
            barcode_ad = adata_subsampled.obs_names[i]
            bc0 = barcode_ad.split("-")[0]
            bc1 = barcode_ad.split("-")[1]
            bc2 = barcode_ad.split("-")[2]
            barcode_infer = bc0+"-"+bc1+"-"+bc2
            g = adata_subsampled.var_names[j]
            infer_subsampled[i][j] = infer.loc[barcode_infer][g]
    return adata_subsampled, infer_subsampled
