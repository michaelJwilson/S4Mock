import numpy as np
import h5py

root  = "/global/project/projectdirs/desi/mocks/bgs/MXXL/small/"
root  = "./"

fpath = root + "galaxy_catalogue_small.hdf5"

print(fpath)

f   = h5py.File(fpath, mode='r')

ra  = f["Data/ra"][...]
dec = f["Data/dec"][...]
z   = f["Data/z_obs"][...]
r   = f["Data/app_mag"][...]
f.close()

print("RA:")
print(ra)

print("Dec:")
print(dec)

print("z:")
print(z)

print("App mag:")
print(r)
