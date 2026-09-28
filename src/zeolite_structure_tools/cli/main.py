import argparse
from pathlib import Path
from ..operations.add_adsorbate import batch_add_adsorbed_molecule

def main():
    p=argparse.ArgumentParser() 
    p.add_argument("--input-dir",type=Path,default=Path("acid_sites"))
    p.add_argument("--output-dir",type=Path,default=Path("output"))
    p.add_argument("--molecule-name",type=str,default="h2o")
    p.add_argument("--Ncore",type=int,default=2)
    a=p.parse_args()
    results = batch_add_adsorbed_molecule(input_dir=a.input_dir,output_dir=a.output_dir,molecule_name=a.molecule_name,Ncore=a.Ncore)
    for r in results:
        if r["success"] == False:
            print(r)

if __name__=="__main__": 
    main()
