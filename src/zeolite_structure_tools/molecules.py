import numpy as np
import pandas as pd
from .geometry import rotate_to_vector,xyz_to_rsf,rotation_back_matrix,rotate

def add_h(direction,distance=0.98):
    v=np.asarray(direction,dtype=float); v/=np.linalg.norm(v); return pd.DataFrame([["H",1,*list(v*distance)]],columns=["type","number","x","y","z"])

def add_mol(direction,template,distance=1.9):
    out=rotate_to_vector(template,direction)
    out[["x","y","z"]]+=np.asarray(direction)*distance
    return out

def normalize_molecule(arc,neighbor_number,neighbor_types=None):
    out=arc.copy()
    origin=out.loc[0,["x","y","z"]].to_numpy(dtype=float)
    out[["x","y","z"]]-=origin
    coords=out[["x","y","z"]].to_numpy(dtype=float)
    out["distance"]=np.linalg.norm(coords,axis=1)
    sel=out if neighbor_types is None else out[out["type"].isin(neighbor_types)]; ids=sel.sort_values("distance").index[:neighbor_number+1]; data=out.loc[ids].copy(); d=data["distance"].to_numpy(); d[d==0]=1; data[["x","y","z"]]=data[["x","y","z"]].to_numpy()/d[:,None]
    if neighbor_number>0:
        v=data[["x","y","z"]].sum().to_numpy(); _,s,f=xyz_to_rsf(v); out[["x","y","z"]]=out[["x","y","z"]].to_numpy()@rotation_back_matrix(-s,-f).T
    return out.drop(columns="distance")
