import os
import sys
from collections import defaultdict
from functools import reduce, partial

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


# data from https://www.bioinformatics.org/sms2/rev_trans.html


def sample():
    N = 7
    data = []
    with open('d.txt') as f:
        for i in range(13):
            f.readline()
            for _ in range(3):
                pos, aa, cdn_pos = f.readline().strip().split('_')
                g, g_share = f.readline().strip().split()
                a, a_share = f.readline().strip().split()
                t, t_share = f.readline().strip().split()
                c, c_share = f.readline().strip().split()
                g_share = float(g_share)
                a_share = float(a_share)
                t_share = float(t_share)
                c_share = float(c_share)

                sample = np.random.choice(
                    list('GATC'), N, 
                    p=[g_share, a_share, t_share, c_share])
                data.append(sample)

    for x in np.array(data).T:
        print(''.join(x)[:18])


color_mapping12 = {
    "C>A": "deepskyblue",
    "G>T": "deepskyblue",
    "C>G": "black",
    "G>C": "black",
    "C>T": "red",
    "G>A": "red",
    "T>A": "silver",
    "A>T": "silver",
    "T>C": "yellowgreen",
    "A>G": "yellowgreen",
    "T>G": "pink",
    "A>C": "pink",
}
sbs12_ordered = ["C>A", "G>T", "C>G", "G>C", "C>T", "G>A", "T>A", "A>T", "T>C", "A>G", "T>G", "A>C"]
colors12 = [color_mapping12[sbs] for sbs in sbs12_ordered]


def plot_mutspec12(
        mutspec: pd.DataFrame, 
        spectra_col="MutSpec", 
        figsize=(4, 8), 
        savepath=None, 
        ticksize=16,
        show=True, 
        **kwargs,
    ):
    fig = plt.figure(figsize=figsize)
    ax = sns.barplot(
        y="Mut", x=spectra_col, data=mutspec, width=0.7, 
        order=sbs12_ordered, ax=fig.gca(), **kwargs)

    # map colors to bars
    for bar, clr in zip(ax.patches, colors12):
        bar.set_color(clr)

    ax.set_title("")
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_xticks([])

    plt.yticks(fontsize=ticksize)

    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()
    return ax


# ["C>A", "G>T", "C>G", "G>C", "C>T", "G>A", "T>A", "A>T", "T>C", "A>G", "T>G", "A>C"]
s = ["4", "3", "6", "5", "30", "17", "8", "6", "10", "20", "4", "3"]
s = np.array([int(x) for x in s])
s = s / s.sum()
print(s)

df = pd.DataFrame({
    'Mut': sbs12_ordered,
    'MutSpec': s,
})

plot_mutspec12(df, savepath='spectrum.pdf', show=False)