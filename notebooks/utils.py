from statistics import geometric_mean

import numpy as np
import pandas as pd
from ete3 import PhyloTree
from pymutspec.constants import possible_sbs192, possible_sbs12
from pymutspec.annotation.tree import get_ingroup_root, get_tree_len, calc_phylocoefs
from pymutspec.annotation.spectra import (
    jackknife_spectra_sampling, calc_edgewise_spectra,
    complete_sbs192_columns, collapse_sbs192, get_cossim, get_eucdist,
)

def assign_cat(p: float, interval=0.1):
    if interval < 0.01:
        raise NotImplementedError
    
    left = p // interval / (1 / interval)
    right = left + interval
    if interval >= 0.1:
        return f"{left:.1f}_{right:.1f}"
    else:
        return f"{left:.2f}_{right:.2f}"


def calc_tree_mutspec(obs: pd.DataFrame, exp: pd.DataFrame, pmin=0.0, pmax=1.0, scale=True):
    obs[obs.ProbaFull.between(pmin, pmax, "right")]
        
    obs = obs[obs.ProbaFull.between(pmin, pmax, "right")]
    exp = exp[exp.Proba.between(pmin, pmax, "right")]
    if obs.ProbaFull.sum() < 100:
        return None
    obs_freqs = obs.groupby(["RefNode", "AltNode", "Mut"]).ProbaFull.sum().unstack()
    assert len(obs_freqs) > 0
    exp_freqs = exp.groupby(["Node", "Mut"]).Proba.sum().unstack()
    assert len(exp_freqs) > 0

    for sbs192 in possible_sbs192:
        if sbs192 not in obs_freqs.columns:
            obs_freqs[sbs192] = 0.0
        if sbs192 not in exp_freqs.columns:
            exp_freqs[sbs192] = 0.0

    edges_df = obs_freqs.index.to_frame(False)
    exp_freqs = edges_df.merge(exp_freqs.reset_index(), left_on="RefNode", right_on="Node")\
        .set_index(["RefNode", "AltNode"])
    
    indexes = obs_freqs.index.intersection(exp_freqs.index)
    obs_freqs = obs_freqs.loc[indexes, possible_sbs192]
    exp_freqs = exp_freqs.loc[indexes, possible_sbs192]
        
    assert (obs_freqs.columns == exp_freqs.columns).all()
    assert (obs_freqs.index == exp_freqs.index).all()

    spectra = (obs_freqs / exp_freqs).fillna(0.)
    if scale:
        spectra = (spectra.T / spectra.T.sum(axis=0)).T

    return spectra
