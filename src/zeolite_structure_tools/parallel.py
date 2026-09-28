from multiprocessing import Pool

def pararun(f, args, Ncore=2):
    with Pool(Ncore) as pool:
        results = pool.starmap(f, args)
    return results