import copy
import pandas as pd

def replace_atoms(arc,indices,element,periodic_table):
    out=arc.copy(); out.loc[list(indices),"type"]=element; out.loc[list(indices),"number"]=periodic_table.loc[element,"number"]; return out

def concat(*frames): return pd.concat([x.copy() for x in frames],ignore_index=True,join="inner")

def deepcopy(arc): return copy.deepcopy(arc)
