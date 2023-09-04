# run from ./data/alisim

import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd

from pymutspec.annotation import CodonAnnotation
from pymutspec.constants import possible_sbs192, possible_sbs12


def complete_sbs192_columns(df: pd.DataFrame):
    df = df.copy()
    if len(df.columns) != 192:
        for sbs192 in set(possible_sbs192).difference(df.columns.values):
            df[sbs192] = 0.
    df = df[possible_sbs192]
    return df


def collapse_sbs192(df: pd.DataFrame, to=12):
    assert (df.columns == possible_sbs192).all()
    df = df.copy()
    if to == 12:
        for sbs192 in possible_sbs192:
            sbs12 = sbs192[2:5]
            if sbs12 in df.columns.values:
                df[sbs12] += df[sbs192]
            else:
                df[sbs12] = df[sbs192]

        return df[possible_sbs12]
    else:
        raise NotImplementedError()


def calc_edgewise_spectra(obs: pd.DataFrame, exp: pd.DataFrame, nmtypes_cutoff=16, nobs_cuttof=20, collapse_to_12=False, scale=True, both_12_and_192=False):
    if len(obs.columns) == 192 and (obs.columns == possible_sbs192).all() and (exp.columns == possible_sbs192).all():
        assert obs.index.names == ["RefNode", "AltNode"]
        assert exp.index.names == ["Node"]
        obs_edges = obs
        freqs_nodes = exp
        freqs_nodes.index.name = "RefNode"
    else:
        obs_edges = obs.groupby(["RefNode", "AltNode", "Mut"]).ProbaFull.sum().unstack()
        obs_edges = complete_sbs192_columns(obs_edges)
        freqs_nodes = exp.groupby(["Node", "Mut"]).Proba.sum().unstack()
        freqs_nodes.index.name = "RefNode"
        freqs_nodes = complete_sbs192_columns(freqs_nodes)

    if not collapse_to_12:
        obs_edges = obs_edges[((obs_edges > 0).sum(axis=1) >= nmtypes_cutoff) & (obs_edges.sum(axis=1) > nobs_cuttof)]
    
    edges_df = obs_edges.index.to_frame(False)

    freqs_edges = edges_df.merge(freqs_nodes, on="RefNode")\
        .set_index(["RefNode", "AltNode"])[possible_sbs192]
    
    
    assert (obs_edges.columns == freqs_edges.columns).all()
    assert (obs_edges.index == freqs_edges.index).all()

    if both_12_and_192:
        obs_edges12   = collapse_sbs192(obs_edges.fillna(0.),   to=12)
        freqs_edges12 = collapse_sbs192(freqs_edges.fillna(0.), to=12)

        spectra12 = (obs_edges12 / freqs_edges12).replace(np.inf, 0.).fillna(0.)
        spectra192 = (obs_edges / freqs_edges).replace(np.inf, 0.).fillna(0.)
        if scale:
            spectra12 = (spectra12.T / spectra12.T.sum(axis=0)).T
            spectra192 = (spectra192.T / spectra192.T.sum(axis=0)).T

        spectra12 = spectra12.fillna(0)
        spectra192 = spectra192.fillna(0)
        assert not (spectra12 == np.inf).any().any()
        assert not (spectra12.isna()).any().any()
        assert not (spectra192 == np.inf).any().any()
        assert not (spectra192.isna()).any().any()
        
        return spectra12, spectra192

    if collapse_to_12:
        obs_edges   = collapse_sbs192(obs_edges.fillna(0.),   to=12)
        freqs_edges = collapse_sbs192(freqs_edges.fillna(0.), to=12)

    spectra = (obs_edges / freqs_edges).replace(np.inf, 0.).fillna(0.)
    if scale:
        spectra = (spectra.T / spectra.T.sum(axis=0)).T

    spectra = spectra.fillna(0)
    assert not (spectra == np.inf).any().any()
    assert not (spectra.isna()).any().any()
    return spectra


# def assign_cat(p: float, interval=0.1):
#     if interval < 0.01:
#         raise NotImplementedError
    
#     left = p // interval / (1 / interval)
#     right = left + interval
#     if interval >= 0.1:
#         return f"{left:.1f}_{right:.1f}"
#     else:
#         return f"{left:.2f}_{right:.2f}"


def get_cossim(a: pd.DataFrame, b: pd.DataFrame):
    assert (a.columns == b.columns).all()
    
    common_index = a.index.intersection(b.index)
    if len(common_index) == 0:
        return pd.Series()
    
    a = a.loc[common_index]
    b = b.loc[common_index]

    dotprod = (a * b).sum(axis=1)
    a_norm = (a ** 2).sum(axis=1) ** 0.5
    b_norm = (b ** 2).sum(axis=1) ** 0.5
    cossim = dotprod / (a_norm * b_norm)
    return cossim


coda = CodonAnnotation(2)

cols_rec_obs = ["Mut", "Label", "PosInGene", "ProbaFull", "RefNode", "AltNode"]
cols_gt_obs  = ["Mut", "Label", "PosInGene", "RefNode", "AltNode"]
cols_rec_exp = ["Mut", "Label", "Pos", "Proba", "Node"]
cols_gt_exp  = ["Mut", "Label", "Pos",  "Node"]


pcutoff = 0.3
cat_cutoff = 4

d = sys.argv[1]
cond = d.split("/")[-1]

path_to_rec_obs = f"spectra_reconstructed_mam/{cond}/spectra_v3/mutations.tsv"
path_to_rec_exp = f"spectra_reconstructed_mam/{cond}/spectra_v3/expected_mutations.tsv"
path_to_gt_obs  = f"spectra_groundtruth_mam/{cond}/mutations.tsv"
path_to_gt_exp  = f"spectra_groundtruth_mam/{cond}/expected_mutations.tsv"
path_to_rates   = os.path.join(d, "IQTREE/anc.rate")

internal_mapping = {f"Node{x}":f"Node{x}" for x in range(5000)}
path_to_mapping = os.path.join(d, "sequences/species_mapping.txt")
if os.path.exists(path_to_mapping):
    name2id_leaves = pd.read_csv(path_to_mapping, header=None, index_col=1, sep="\t")[0].to_dict()
    name2id = dict(**internal_mapping, **name2id_leaves)
else:
    name2id = internal_mapping

if os.path.exists(path_to_rec_exp) and os.path.exists(path_to_gt_exp):
    rec_obs = pd.read_csv(path_to_rec_obs, sep="\t", usecols=cols_rec_obs)
    rec_exp = pd.read_csv(path_to_rec_exp, sep="\t", usecols=cols_rec_exp)
    gt_obs  = pd.read_csv(path_to_gt_obs,  sep="\t", usecols=cols_gt_obs)
    gt_exp  = pd.read_csv(path_to_gt_exp,  sep="\t", usecols=cols_gt_exp)
    rates   = pd.read_csv(path_to_rates,   sep="\t", comment="#")

    rec_obs = rec_obs.merge(rates[["Site", "Cat"]], left_on="PosInGene", right_on="Site").drop("PosInGene", axis=1)
    rec_exp = rec_exp.merge(rates[["Site", "Cat"]], left_on="Pos", right_on="Site").drop("Pos", axis=1)
    gt_obs  = gt_obs.merge(rates[["Site", "Cat"]], left_on="PosInGene", right_on="Site").drop("PosInGene", axis=1)
    gt_exp  = gt_exp.merge(rates[["Site", "Cat"]], left_on="Pos", right_on="Site").drop("Pos", axis=1)

    # replace encoded node names
    gt_obs["AltNode"] = gt_obs.AltNode.map(name2id)

    #prepare OBS
    rec_obs_all = rec_obs[(rec_obs.ProbaFull > pcutoff) & (rec_obs.Label >= 0) & (rec_obs.Cat >= cat_cutoff)]
    rec_obs_syn = rec_obs[(rec_obs.ProbaFull > pcutoff) & (rec_obs.Label >= 1) & (rec_obs.Cat >= cat_cutoff)]
    gt_obs_all = gt_obs[(gt_obs.Label >= 0) & (gt_obs.Cat >= cat_cutoff)].assign(ProbaFull=1.0)
    gt_obs_syn = gt_obs[(gt_obs.Label >= 1) & (gt_obs.Cat >= cat_cutoff)].assign(ProbaFull=1.0)

    #prepare EXP
    rec_exp_all = rec_exp[(rec_exp.Proba > pcutoff) & (rec_exp.Label.astype("category") == "all") & (rec_exp.Cat >= cat_cutoff)]
    rec_exp_syn = rec_exp[(rec_exp.Proba > pcutoff) & (rec_exp.Label.astype("category") == "syn") & (rec_exp.Cat >= cat_cutoff)]
    gt_exp_all = gt_exp[(gt_exp.Label.astype("category") == "all") & (gt_exp.Cat >= cat_cutoff)].assign(Proba=1.0)
    gt_exp_syn = gt_exp[(gt_exp.Label.astype("category") == "syn") & (gt_exp.Cat >= cat_cutoff)].assign(Proba=1.0)
    
    #calc tree spectra
    rec_all_spectra12, rec_all_spectra192 = calc_edgewise_spectra(rec_obs_all, rec_exp_all, 0, 0, both_12_and_192=True)
    gt_all_spectra12, gt_all_spectra192   = calc_edgewise_spectra(gt_obs_all, gt_exp_all, 0, 0, both_12_and_192=True)
    rec_syn_spectra12, rec_syn_spectra192 = calc_edgewise_spectra(rec_obs_syn, rec_exp_syn, 0, 0, both_12_and_192=True)
    gt_syn_spectra12, gt_syn_spectra192   = calc_edgewise_spectra(gt_obs_syn, gt_exp_syn, 0, 0, both_12_and_192=True)
    
    #get cossim
    cossim12_all  = get_cossim(gt_all_spectra12,  rec_all_spectra12)
    cossim192_all = get_cossim(gt_all_spectra192, rec_all_spectra192)
    cossim12_syn  = get_cossim(gt_syn_spectra12,  rec_syn_spectra12)
    cossim192_syn = get_cossim(gt_syn_spectra192, rec_syn_spectra192)

    total_cossim = cossim12_all.merge(cossim12_syn, "outer", on=["RefNode", "AltNode"])\
        .merge(cossim192_all, "outer", on=["RefNode", "AltNode"])\
            .merge(cossim192_syn, "outer", on=["RefNode", "AltNode"])
    
    total_cossim.to_csv("../")