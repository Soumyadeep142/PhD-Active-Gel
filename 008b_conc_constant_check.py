from numpy import *
import glob
import os

files = sorted(glob.glob("step_*.dat"), key=lambda x: int(x.split("_")[1].split(".")[0]))
for f in files:
    data = loadtxt(f)
    x, y, z, phi, conc = data[:,0], data[:,1], data[:,2], data[:,3],data[:,4]
    
    print(f, sum(conc))
