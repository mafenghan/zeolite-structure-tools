from pathlib import Path
import numpy as np
from ..io import read_arc
from ..geometry import distance_matrix_pbc, find_neighbors, abc2latt, car2dir, dir2car, find_non_collision_rotation
from ..molecules import normalize_molecule, add_mol
from ..structure import concat
from ..arc_writer import write_arc
from ..parallel import pararun

def read_template_molecule(molecule_name):
    if molecule_name.lower() == "h2o":
        template_path = Path('data/templates/H2O.arc')
        cn = 2
    elif molecule_name.lower() == "co":
        template_path = Path('data/templates/CO.arc')
        cn = 1
    elif molecule_name.lower() == "c2h4o":
        template_path = Path('data/templates/C2H4O.arc')
        cn = 2
    elif molecule_name.lower() == "nh3":
        template_path = Path('data/templates/NH3.arc')
        cn = 3
    elif molecule_name.lower() == "c2h5n":
        template_path = Path('data/templates/C2H5N.arc')
        cn = 3
    elif molecule_name.lower() == "c2h4":
        template_path = Path('data/templates/C2H4.arc')
        cn = 3
    elif molecule_name.lower() == "c3h6":
        template_path = Path('data/templates/C3H6.arc')
        cn = 4
    else:
        raise ValueError(f"Unsupported molecule: {molecule_name}")
    template_path = Path(template_path)
    if not template_path.exists():
        raise FileNotFoundError(f"Template file not found: {template_path}")
    return normalize_molecule(read_arc(template_path).arc, cn)

def add_adsorbed_molecule(name,input_dir,output_dir="output",molecule_name="h2o",collision_threshold=1.4,placement_length=2.7):
    input_dir=Path(input_dir)
    output_dir=Path(output_dir)
    output_dir.mkdir(parents=True,exist_ok=True)

    template = read_template_molecule(molecule_name)
    st=read_arc(input_dir/f"{name}.arc")
    # write_arc(template, st.abc, output_dir/f"{molecule_name.lower()}_template.arc", sort=False)
    st.lattice = abc2latt(st.abc)
    arc=st.arc
    arc[["x","y","z"]]=car2dir(st.lattice,arc[["x","y","z"]].to_numpy(dtype=float))
    T=arc[arc["type"].isin(["Si","Al","P"])].reset_index(drop=True)
    O=arc[arc["type"]=="O"].reset_index(drop=True)
    H=arc[arc["type"]=="H"].reset_index(drop=True)
    if H.empty: return {"host":name,'guest':molecule_name,"success":False,"reason":"no H"}
    d_oh=distance_matrix_pbc(O,H,st.lattice)
    oi=np.where(d_oh[:,0]<1.05)[0]
    if len(oi)==0: return {"host":name,'guest':molecule_name,"success":False,"reason":"no O-H site"}
    oi=oi[0]
    d_to = distance_matrix_pbc(T,O,st.lattice)
    ti = find_neighbors(d_to,oi,2.5)
    if len(ti)==0: return {"host":name,'guest':molecule_name,"success":False,"reason":"no T-O site"}
    v1=np.asarray(H.loc[0,["x","y","z"]]-O.loc[oi,["x","y","z"]],float)
    v1-=np.round(v1)
    v1/=np.linalg.norm(v1)
    v2=np.asarray(T.loc[ti[0],["x","y","z"]]-O.loc[oi,["x","y","z"]],float)
    v2-=np.round(v2)
    vn=np.cross(v1,v2)
    vn/=np.linalg.norm(vn)
    for m in (0,-0.1,0.1,0.2,-0.2,0.3,-0.3,0.4,-0.4):
        v3=v1+m*vn
        v3 = dir2car(st.lattice, v3.reshape(1,3)).flatten()
        v3/=np.linalg.norm(v3)
        mol=add_mol(v3,template,distance=placement_length)
        rotation_angles = range(0, 271, 30)
        result = find_non_collision_rotation(
            mol=mol,
            environment_arc=arc,
            lattice=st.lattice,
            rotation_axis=v3,
            rotation_angles=rotation_angles,
            anchor_position=O.loc[oi, ["x", "y", "z"]].to_numpy(dtype=float),
            collision_threshold=collision_threshold,
        )
        if result is not None:
            out = concat(arc, result)
            out[["x","y","z"]]=dir2car(st.lattice,out[["x","y","z"]].to_numpy(dtype=float))
            idH= out[out["type"]=="H"].index[0]
            idX = arc.shape[0]
            idF= out[out["type"]=="F"].index
            out.loc[idF,"type"]="H"
            output=output_dir/f"{name}_add{molecule_name.lower()}_{idH}_{idX}.arc"
            write_arc(out,st.abc,output,sort=False)
            return {"host":name,'guest':molecule_name,"success":True,"output":str(output)}
    return {"host":name,'guest':molecule_name,"success":False,"reason":"no collision-free placement"}

def batch_add_adsorbed_molecule(input_dir,output_dir="output",molecule_name="h2o",collision_threshold=1.4,placement_length=2.7, Ncore=2):
    input_dir=Path(input_dir)
    output_dir = Path(output_dir)
    data = [(f.stem,input_dir,output_dir,molecule_name,collision_threshold,placement_length) for f in input_dir.glob("*.arc")]
    results = pararun(add_adsorbed_molecule, data, Ncore=Ncore)
    return results


def add_co(input_dir,output_dir="output",molecule_name="co"):
    input_dir=Path(input_dir)
    output_dir=Path(output_dir)
    output_dir.mkdir(parents=True,exist_ok=True)

    template = read_template_molecule(molecule_name)
    st=read_arc(input_dir/f"MnO.arc")
    st.lattice = abc2latt(st.abc)
    arc=st.arc
    arc[["x","y","z"]]=car2dir(st.lattice,arc[["x","y","z"]].to_numpy(dtype=float))
    T=arc[arc["type"].isin(["Mn","Al","P"])].reset_index(drop=True)
    O=arc[arc["type"]=="O"].reset_index(drop=True)
    d_oT=distance_matrix_pbc(O,T,st.lattice)
    for i in T.index:
        ti = find_neighbors(d_oT,i,2.4)
        print(f"Mn{i}: cn = {len(ti)}")
        if len(ti) > 0:
            v1=np.asarray(O.loc[ti,["x","y","z"]]-T.loc[i,["x","y","z"]],float)
            v1 = -np.sum(v1,axis=0)
            v1-=np.round(v1)
            v1 = dir2car(st.lattice, v1.reshape(1,3)).flatten()
            v1/=np.linalg.norm(v1)
            mol=add_mol(v1,template,distance=1.8)
            result = find_non_collision_rotation(
                mol=mol,
                environment_arc=arc,
                lattice=st.lattice,
                rotation_axis=v1,
                rotation_angles=[0],
                anchor_position=T.loc[i, ["x", "y", "z"]].to_numpy(dtype=float),
                collision_threshold=1.6,
            )
            if result is not None:
                out = concat(arc, result)
                out[["x","y","z"]]=dir2car(st.lattice,out[["x","y","z"]].to_numpy(dtype=float))
                output=output_dir/f"MnO_add{molecule_name.lower()}_{i}.arc"
                write_arc(out,st.abc,output,sort=False)


if __name__=="__main__":
    # data = add_adsorbed_molecule("SAPO_AWO_1T4O14",input_dir="acid_sites", molecule_name = 'C2H5N')
    # print(data)
    add_co(input_dir="acid_sites", molecule_name = 'CO')

