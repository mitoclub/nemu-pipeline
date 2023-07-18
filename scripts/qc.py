import numpy as np
import pandas as pd

DIR = "~/nemu-pipeline/data/bacteria/work/d2/182bbfea5558ef142bdb6f00b1da37/"

leaves = pd.read_csv(DIR + "leaves_states.state", sep='\t')
anc = pd.read_csv(DIR + "iqtree_anc.state", sep='\t')

print(leaves)