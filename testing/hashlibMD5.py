#!/usr/bin/env python
import hashlib
import netCDF4 as nc4

def update_hash_var(var_obj, hash_obj):
    """Update `hash_obj` from a variable."""
    hash_obj.update(str(var_obj.name).encode('utf8'))
    hash_obj.update(var_obj[:].data)
    hash_obj.update(var_obj[:].data)
    return hash_obj

def calculate_file_hash(file_name):
    """Calculate hash for a given netCDF file."""
    hashList=[]
    with nc4.Dataset(str(file_name)) as data_set:
        hash_obj = hashlib.md5()
        for key, var in sorted(data_set.variables.items()):
            hash_obj = update_hash_var(var, hash_obj)
            hashList.append(hash_obj.hexdigest())
        return hashList

if __name__ == "__main__":
    import argparse
    #import glob
    parser = argparse.ArgumentParser(description="generate md5 hash of variables in netcdf file")
    parser.add_argument("-f",type=str,required=True,help="name of netcdf file",nargs="+")
    args = parser.parse_args()
    #for name in glob.glob(args.f):
    for name in args.f:
        print(calculate_file_hash(name), name)
