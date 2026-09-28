from pathlib import Path
import pandas as pd
import re
from .config import Config

class arcread():
    def __init__(self,filelist,sort = False):
        cfg = Config('config/default.yaml')
        PT = read_periodic_table(cfg.periodic_table)
        self.arc = []
        for i in filelist:
            if 'CORE' in i or 'MOL' in i:
                a = i.split()
                atomtype = re.split(r'\d+', a[0])[0]
                if atomtype == 'xx' or '_' in atomtype:
                    atomtype = a[-2]
                self.arc.append([atomtype,PT.loc[atomtype,'number'],float(a[1]),float(a[2]),float(a[3]),int(a[5])])
            elif 'PBC  ' in i:
                a = i.split()
                self.abc = [float(a[1]),float(a[2]),float(a[3]),float(a[4]),float(a[5]),float(a[6])]
            elif 'Energy' in i:
                a = i.split()
                try:
                    self.en =  float(a[3])
                except:
                    self.en =  0.0
        self.arc = pd.DataFrame(self.arc,columns=['type','number','x','y','z','charge'])
        if 0 not in self.arc['charge']:
            self.arc['charge'] = [0] * self.arc.shape[0]
        if sort == True:
            self.arc.sort_values(by=['number','z'],inplace = True,ascending = [False,True],ignore_index = True)
    def atominfo(self):
        self.atomtype, self.atomnumber = [], []
        for i in self.arc.index:
            if self.arc.loc[i,'type'] not in self.atomtype:
                self.atomtype.append(self.arc.loc[i,'type'])
                self.atomnumber.append(1)
            else:
                self.atomnumber[-1] += 1

def read_arc(path_or_lines):
    if isinstance(path_or_lines,(str,Path)):
        path_or_lines=Path(path_or_lines).read_text(encoding="utf-8").splitlines()
    return arcread(path_or_lines)

def read_periodic_table(path):
    return pd.read_csv(Path(path),index_col=0)
