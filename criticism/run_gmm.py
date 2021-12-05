from cmdstanpy import model
import numpy as np
import scanpy as sc
import json
import pandas as pd
import os
from cmdstanpy import cmdstan_path, CmdStanModel
import scipy.stats as stats
import matplotlib.pyplot as plt
import random

from cmdstanpy import cmdstan_path, set_cmdstan_path
set_cmdstan_path('/Users/alexander/cmdstan')


PATH_TO_STAN_FILES = "/Users/alexander/classes/Ml prob prog/project/"




#Create the data object
def create_gmm_data(small_infer):
    data = {}
    data['C'] = len(small_infer)
    data['G'] = len(small_infer[0])


    data['prior_mu'] = [.9109,1.001,1.1167]

    data['infer'] = small_infer

    GMM_data = os.path.join(cmdstan_path(), 'scHPF_GMM','scHPF_GMM.data.json')

    with open(GMM_data, "w") as outfile:
        json.dump(data, outfile)
    return GMM_data






def run_gmm(stan_model,infer):
    GMM_data = create_gmm_data(infer)
    print(PATH_TO_STAN_FILES)
    print(stan_model)
    GMM_stan = os.path.join(PATH_TO_STAN_FILES,stan_model)
    print(GMM_stan)
    GMM_model = CmdStanModel(stan_file=GMM_stan)
    GMM_model.name
    GMM_model.stan_file
    GMM_model.exe_file
    GMM_model.code()
    GMM_fit = GMM_model.variational(data=GMM_data, output_dir='.',show_console=True,require_converged=False,tol_rel_obj=.001,iter=20000)
    return GMM_fit
