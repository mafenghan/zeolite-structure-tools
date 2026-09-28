from pathlib import Path

def write_arc(arc,abc,filename,sort=True):
    data=arc.copy()
    if sort:
        data=data.sort_values(by=["number","z"],ascending=[False,True])
    out=Path(filename); out.parent.mkdir(parents=True,exist_ok=True)
    with out.open("w",encoding="utf-8") as f:
        f.write("!BIOSYM archive 2\nPBC=ON\n         Energy   0    0      0.000000\n!DATE\n")
        f.write("PBC   %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n"%tuple(abc[:6]))
        for n,(_,r) in enumerate(data.iterrows(),1):
            t=str(r["type"]); fmt=("%s %18.9f %14.9f %14.9f %s %4s %2s %2s %4s %4s\n" if len(t)==1 else "%s %17.9f %14.9f %14.9f %s %4s %2s %2s %4s %4s\n")
            f.write(fmt%(t,r["x"],r["y"],r["z"],"CORE",n,t,t,"0.0000",n))
        f.write("end\nend\n")
