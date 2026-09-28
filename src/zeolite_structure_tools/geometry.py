import numpy as np
from math import pi, sqrt, sin, cos

def latt2abc(latt):
    abc = np.zeros((6))
    abc[0] = np.linalg.norm(latt[0,:])
    abc[1] = np.linalg.norm(latt[1,:])
    abc[2] = np.linalg.norm(latt[2,:])
    ua = (latt[0,:]/abc[0])
    ub = (latt[1,:]/abc[1])
    uc = (latt[2,:]/abc[2])
    abc[3] = 180.0/pi*np.arccos(np.dot(ub,uc))
    abc[4] = 180.0/pi*np.arccos(np.dot(ua,uc))
    abc[5] = 180.0/pi*np.arccos(np.dot(ua,ub))
    return abc

def abc2latt(abc):
    a = abc[0]
    b = abc[1]
    c = abc[2]
    alf = abc[3]
    bet = abc[4]
    gam = abc[5]
    abc_latt = np.zeros((3,3))
    abc_latt[0,0] = a
    abc_latt[1,0] = b*cos(gam*pi/180.0)
    abc_latt[1,1] = b*sin(gam*pi/180.0)
    abc_latt[2,0] = c*cos(bet*pi/180.0)
    abc_latt[2,1] = (b*c*cos(alf*pi/180.0)-abc_latt[1,0]*abc_latt[2,0])/abc_latt[1,1]
    abc_latt[2,2] = sqrt(c**2-abc_latt[2,0]**2-abc_latt[2,1]**2)
    return np.mat(abc_latt)

def dir2car(lattice,atomdir):
    atomcar = np.zeros((atomdir.shape[0],3))
    for i in range(atomdir.shape[0]):
        atomcar[i] = np.mat(atomdir[i,:])*np.mat(lattice)
    return np.array(atomcar)

def car2dir(lattice,atomcar):
    atomdir = np.mat(atomcar) * np.mat(lattice).I
    return np.array(atomdir)

def distance_matrix(row,col):
    a=row[["x","y","z"]].to_numpy(dtype=float); b=col[["x","y","z"]].to_numpy(dtype=float)
    return np.linalg.norm(b[:,None,:]-a[None,:,:],axis=2)

def distance_matrix_pbc(r,c,lattice):
    r_max_xyz = r[['x','y','z']].to_numpy(dtype=float).max()
    c_max_xyz = c[['x','y','z']].to_numpy(dtype=float).max()
    max_xyz = max(r_max_xyz, c_max_xyz)
    if max_xyz > 1.0:
        r_frac_coords = car2dir(lattice, r[['x','y','z']].to_numpy(dtype=float))
        c_frac_coords = car2dir(lattice, c[['x','y','z']].to_numpy(dtype=float))
    else:
        r_frac_coords = r[['x','y','z']].to_numpy(dtype=float)
        c_frac_coords = c[['x','y','z']].to_numpy(dtype=float)
    r_frac_coords = np.mod(r_frac_coords, 1.0)
    c_frac_coords = np.mod(c_frac_coords, 1.0)
    dist = []
    for i in range(c_frac_coords.shape[0]):
        delta = np.array(r_frac_coords - c_frac_coords[i,:],dtype = float)
        delta -= np.round(delta) # 去除周期性的影响
        dist.append(np.linalg.norm(np.dot(delta, lattice),axis = 1))
    return np.array(dist).T

def find_neighbors(distances,atom_index,cutoff,min_distance=0.5):
    v=distances[:,atom_index]; return [i for i,x in enumerate(v) if min_distance<x<cutoff]

def xyz_to_rsf(v):
    v=np.asarray(v,dtype=float); r=np.linalg.norm(v)
    if r==0: raise ValueError("zero vector has no direction")
    x,y,z=v/r; return np.array([r,np.arctan2(y,x),np.arctan2(np.hypot(x,y),z)])

def rotation_matrix(s,f):
    Ry=np.array([[np.cos(f),0,np.sin(f)],[0,1,0],[-np.sin(f),0,np.cos(f)]])
    Rz=np.array([[np.cos(s),-np.sin(s),0],[np.sin(s),np.cos(s),0],[0,0,1]])
    return Rz@Ry

def rotation_back_matrix(s,f):
    Ry=np.array([[np.cos(f),0,np.sin(f)],[0,1,0],[-np.sin(f),0,np.cos(f)]])
    Rz=np.array([[np.cos(s),-np.sin(s),0],[np.sin(s),np.cos(s),0],[0,0,1]])
    return Ry@Rz

def rotate(arc,R):
    out=arc.copy(); out[["x","y","z"]]=arc[["x","y","z"]].to_numpy(dtype=float)@R.T; return out

def rotate_to_vector(arc,v):
    r,s,f=xyz_to_rsf(v); return rotate(arc,rotation_matrix(s,f))

def rotation_matrix_axis(v,angle_deg):
    v=np.asarray(v,dtype=float); v=v/np.linalg.norm(v); a=np.deg2rad(angle_deg); c=np.cos(a); s=np.sin(a); x,y,z=v
    return np.array([[x*x*(1-c)+c,x*y*(1-c)-z*s,x*z*(1-c)+y*s],[x*y*(1-c)+z*s,y*y*(1-c)+c,y*z*(1-c)-x*s],[x*z*(1-c)-y*s,y*z*(1-c)+x*s,z*z*(1-c)+c]])

def rotate_axis(arc,v,angle_deg): return rotate(arc,rotation_matrix_axis(v,angle_deg))

def supercell(arc,abc):
    lattice=abc2latt(abc); atomdir=car2dir(lattice,arc[["x","y","z"]].to_numpy(dtype=float)); rows=[]
    for i in range(len(atomdir)):
        for j in (-1,0,1):
            for k in (-1,0,1):
                for q in (-1,0,1): rows.append([arc.loc[i,"type"],i,atomdir[i,0]+j,atomdir[i,1]+k,atomdir[i,2]+q])
    import pandas as pd
    out=pd.DataFrame(rows,columns=["type","relation","x","y","z"])
    if min(abc[:3])>4: out=out[(out[["x","y","z"]].max(axis=1)<=1.5)&(out[["x","y","z"]].min(axis=1)>=-0.5)].reset_index(drop=True)
    out[["x","y","z"]]=dir2car(lattice,out[["x","y","z"]].to_numpy()); return out

def find_non_collision_rotation(
    mol,
    environment_arc,
    lattice,
    rotation_axis,
    rotation_angles,
    anchor_position,
    collision_threshold,
):
    """
    旋转分子并寻找一个与周围原子不发生碰撞的构型。

    Parameters
    ----------
    mol : pd.DataFrame
        待添加的分子，坐标为笛卡尔坐标或当前函数要求的坐标体系。

    environment_arc : pd.DataFrame
        环境的原子结构，坐标为分数坐标。

    lattice
        晶格信息。

    rotation_axis
        分子旋转轴。

    rotation_angles : iterable
        待尝试的旋转角度，例如 range(0, 271, 30)。

    anchor_position
        分子锚定位置，坐标为分数坐标。

    collision_threshold : float
        碰撞判断距离阈值。

    Returns
    -------
    tuple or None
        成功：
            (rotated_mol, angle)

        如果所有角度都发生碰撞：
            None
    """
    anchor_position = np.asarray(anchor_position, dtype=float)
    for angle in rotation_angles:
        # 1. 旋转分子
        rotated = rotate_axis(
            mol,
            rotation_axis,
            angle,
        )
        # 2. 转换到 fractional/direct 坐标
        rotated[["x", "y", "z"]] = car2dir(
            lattice,
            rotated[["x", "y", "z"]].to_numpy(dtype=float),
        )
        # 3. 将分子平移到锚定位置
        rotated[["x", "y", "z"]] += anchor_position
        # 4. 判断是否发生碰撞
        dm = distance_matrix_pbc(
            environment_arc,
            rotated,
            lattice,
        )
        if not np.any(dm < collision_threshold):
            return rotated
    # 所有角度均发生碰撞
    return None


