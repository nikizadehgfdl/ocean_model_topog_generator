#!/usr/bin/env python
import GMesh_torch
import netCDF4
import numpy as np
import torch
import argparse
import time

def convol( levels, h, f, verbose=False ):
    """Coarsens the product of h*f across all levels"""
    levels[-1].height = ( h * f ).reshape(levels[-1].nj,levels[-1].ni)
    for k in range( len(levels) - 1, 0, -1 ):
        if verbose: print('Coarsening {} -> {}'.format(k,k-1))
        levels[k].coarsenby2( levels[k-1] )
    return levels[0].height  

def rough( levels, h, h2min=1.e-7 ):
    """Calculates both mean of H, and variance of H relative to a plane"""
    # Construct weights for moment calculations
    nx = 2**( len(levels) - 1 )
    x = ( np.arange(nx) - ( nx - 1 ) /2 ) * np.sqrt( 12 / ( nx**2 - 1 ) ) # This formula satisfies <x>=0 and <x^2>=1
    X, Y = np.meshgrid( x, x )
    X, Y = X.reshape(1,nx,1,nx), Y.reshape(1,nx,1,nx)
    h = h.reshape(levels[0].nj,nx,levels[0].ni,nx)
    # Now calculate moments
    H2 = convol( levels, h, h ) # mean of h^2
    HX = convol( levels, h, X ) # mean of h * x
    HY = convol( levels, h, Y ) # mean of h * y
    H = convol( levels, h, np.ones((1,nx,1,nx)) ) # mean of h = mean of h * 1
    # The variance of deviations from the plane = <h^2> - <h>^2 - <h*x>^2 - <h*y>^2 given <x>=<y>=0 and <x^2>=<y^2>=1
    return H, H2 - H**2 - HX**2 - HY**2 + h2min

def rough_opt1( levels, h , device, h2min=1.e-7):
    """Calculates both mean of H, and variance of H relative to a plane"""
    # Construct weights for moment calculations
    nx = 2**( len(levels) - 1 )
    x = ( np.arange(nx) - ( nx - 1 ) /2 ) * np.sqrt( 12 / ( nx**2 - 1 ) ) # This formula satisfies <x>=0 and <x^2>=1
    X=torch.from_numpy(X).to(device)
    Y=torch.from_numpy(Y).to(device)
    h = torch.reshape(h,(levels[0].nj,nx,levels[0].ni,nx)).to(device) #Why is this not already on device? Input h is on device!
    # Now calculate moments
    H2 = convol( levels, h, h ) # mean of h^2
    HX = convol( levels, h, X ) # mean of h * x
    HY = convol( levels, h, Y ) # mean of h * y
    H = convol( levels, h, torch.ones((1,nx,1,nx)).to(device)) # mean of h = mean of h * 1
    # The variance of deviations from the plane = <h^2> - <h>^2 - <h*x>^2 - <h*y>^2 given <x>=<y>=0 and <x^2>=<y^2>=1
    return H, H2 - H**2 - HX**2 - HY**2 + torch.tensor((h2min))

def do_RSC_new(targG,src_topo_global, NtileI=1, NtileJ=1, max_refinement=10, device='cpu'):
    """Apply the RSC algoritm using a fixed number of refinements max_refinement"""
    di, dj = targG.ni // NtileI, targG.nj // NtileJ
    assert di*NtileI == targG.ni
    assert dj*NtileJ == targG.nj
    print('window size dj,di =',dj,di,'full model nj,ni=',targG.nj, targG.ni)
    Hcnt = np.zeros((targG.nj, targG.ni)) # Diagnostic: counting which cells we are working on
    Htarg, H2targ = np.zeros((targG.nj, targG.ni)), np.zeros((targG.nj, targG.ni))
    for j in range(NtileJ ):
        csj, sj = slice( j*dj, (j+1)*dj ), slice( j*dj, (j+1)*dj+1 )
        for i in range(NtileI ):
            csi, si = slice( i*di, (i+1)*di ), slice( i*di, (i+1)*di+1 ) # Slices of target grid
            Hcnt[csj,csi] = Hcnt[csj,csi] + 1 # Diagnostic: counting which cells we are working on
            G = GMesh_torch.GMesh_torch( lon=targG.lon[sj,si], lat=targG.lat[sj,si], device=device)
            print('J,I={},{} {:.1f}%, {}\n   window lon={}:{}, lat={}:{}\n   jslice={}, islice={}'.format( \
                j, i, 100*(j*NtileI+i)/(NtileI*NtileJ), G, G.lon.min(), G.lon.max(), G.lat.min(), G.lat.max(), sj, si ))
            levels = G.refine_loop(src_topo_global, resolution_limit=False, fixed_refine_level=max_refinement, 
                                   timers=False, verbose=False, device=device)
            ## Use nearest neighbor topography to populate the finest grid
            levels[-1].project_source_data_onto_target_mesh( src_topo_global)
            h, h2 = rough(levels, levels[-1].height)
            # Store window in final array
            Htarg[csj,csi] = h
            H2targ[csj,csi] = h2
    print( Hcnt.min(), Hcnt.max(), '<-- should both be 1 for full model' )
    return  Htarg, H2targ
    
def write_topog(targH,targH2,targlon,targlat,filename,description=None,history=None,source=None,no_changing_meta=None):
    if filename is None:
        filename='topog.nc'
    with netCDF4.Dataset(filename,'w','clobber') as nc:
        nx = nc.createDimension('nx', targH.shape[1])
        ny = nc.createDimension('ny', targH.shape[0])
        ntiles = nc.createDimension('ntiles', 1)
        z = nc.createVariable('depth', float, ('ny','nx'))
        z.units='meters' 
        z2 = nc.createVariable('h2', float, ('ny','nx'))
        z2.units='meters^2'
        z[:,:] = -targH[:,:]
        z2[:,:] = targH2[:,:]
        #global attributes
        if(not no_changing_meta):
            nc.history = history
            nc.description = description
            nc.source =  source

        #x=nc.createVariable('x','f8',('ny','nx'))
        #x.units='meters'
        #x[:]=xx
        #y=nc.createVariable('y','f8',('ny','nx'))
        #y.units='meters'
        #y[:]=yy

def global_meta_info():
    import socket, os, subprocess, datetime, sys
    
    host = str(socket.gethostname())
    scriptpath = sys.argv[0]
    scriptbasename = (subprocess.check_output("basename " + scriptpath, shell=True).decode("ascii").rstrip("\n"))
    scriptdirname = (
        subprocess.check_output("dirname " + scriptpath, shell=True)
        .decode("ascii")
        .rstrip("\n")
    )
    scriptgithash = (
        subprocess.check_output(
            "cd " + scriptdirname + ";git rev-parse HEAD; exit 0",
            stderr=subprocess.STDOUT,
            shell=True,
        )
        .decode("ascii")
        .rstrip("\n")
    )
    scriptgitMod = (
        subprocess.check_output(
            "cd "
            + scriptdirname
            + ";git status --porcelain "
            + scriptbasename
            + " | awk '{print $1}' ; exit 0",
            stderr=subprocess.STDOUT,
            shell=True,
        )
        .decode("ascii")
        .rstrip("\n")
    )
    if "M" in str(scriptgitMod):
        scriptgitMod = " , But was localy Modified!"

    hist = "This file was generated via command " + " ".join(sys.argv)
    hist = hist + " on " + str(datetime.date.today()) + " on platform " + host

    desc = ("This file contains the ocean depth for the input model grid. ")

    source = ""
    source = source + scriptpath + " had git hash " + scriptgithash + scriptgitMod
    source = (source
        + ". To obtain the grid generating code do: git clone  https://github.com/nikizadehgfdl/ocean_model_topog_generator.git ; cd ocean_model_topog_generator;  git checkout "
        + scriptgithash
    )

    return hist,desc,source
        
def main(hgridfilename,outputfilename,
         source_file,source_lon,source_lat,source_elv,
         device,nxblocks,nyblocks,max_refine,no_changing_meta):

    desc=''
    hist=''
    source=''
    host=''
    if not no_changing_meta:
        hist,desc,source=global_meta_info()

    #Time it
    tic = time.perf_counter()
    #Open and read the topographic dataset
    with netCDF4.Dataset(source_file) as nc:
        topo_lon = nc.variables[source_lon][:].filled(0.)
        topo_lat = nc.variables[source_lat][:].filled(0.)
        topo_elv = nc.variables[source_elv][:,:].filled(0.)
    
    #Create the topo object that contains the source data
    t_topo_lon=torch.from_numpy(topo_lon).to(device)
    t_topo_lat=torch.from_numpy(topo_lat).to(device)
    t_topo_elv=torch.from_numpy(topo_elv).to(device)
    src_topo_global = GMesh_torch.UniformEDS( t_topo_lon, t_topo_lat, t_topo_elv, device)
    print("source shape: ",src_topo_global.data.shape)
    print("source lon range: ",src_topo_global.lonh[0],src_topo_global.lonh[-1])
    print("source lat range: ",src_topo_global.lath[0],src_topo_global.lath[-1])
    
    #Open and read a target grid
    with netCDF4.Dataset(hgridfilename) as nc:
        lon=nc.variables['x'][::2,::2]
        lat=nc.variables['y'][::2,::2]
    #create the target mesh object that contains the target grid
    t_lon=torch.from_numpy(lon).to(device)
    t_lat=torch.from_numpy(lat).to(device)
    targG = GMesh_torch.GMesh_torch( lon=t_lon, lat=t_lat, device=device)

    print("target shape: ",targG.lon.shape)
    print("target lon range: ",targG.lon[0,0],targG.lon[0,-1])
    print("target lat range: ",targG.lat[0,0],targG.lat[-1,0])

    #Do the RSC algorithm to deduce depth and roughness
    st = time.time()
    Htarg, H2targ = do_RSC_new(targG,src_topo_global, nxblocks, nyblocks, max_refine, device=device)    
    et = time.time() - st
    print('RSC loop time:', et, 'seconds')
    
    #Write the ocean_topog.nc file
    if (outputfilename is None ):
        outputfilename='new_topo_OM5_grid_r{}_{}x{}.nc'.format(max_refine,nxblocks, nyblocks )
    write_topog(Htarg,H2targ,targG.lon,targG.lat,filename=outputfilename,description=desc,history=hist,source=source,no_changing_meta=no_changing_meta)

    toc = time.perf_counter()
    print(f"It took {toc - tic:0.4f} seconds on platform ",host)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="create ocean topog")

    parser.add_argument("--hgridfilename",type=str,required=False,default="ocean_hgrid.nc",
                        help="pathname of input hgrid file")
    parser.add_argument("--outputfilename",type=str,required=False,default=None,
                        help="pathname of output topog file")
    parser.add_argument("--source_file",type=str,required=False,default="GEBCO_2023.nc",
                        help="pathname of source topography data file")
    parser.add_argument("--source_lon",type=str,required=False,default="lon",
                        help="pathname of longitude variable in source topography data file")
    parser.add_argument("--source_lat",type=str,required=False,default="lat",
                        help="pathname of latitude variable in source topography data file")
    parser.add_argument("--source_elv",type=str,required=False,default="elevation",
                        help="pathname of elevetion variable in source topography data file")
    parser.add_argument("--device",type=str,required=False,default="cpu",
                        help="Specify if gpu should be used by passing --device cuda")
    parser.add_argument("--nxblocks",type=int,required=False,default=1,
                        help="number of i-direction blocks to be used")
    parser.add_argument("--nyblocks",type=int,required=False,default=1,
                        help="number of j-direction blocks to be used")
    parser.add_argument("--max_refine",type=int,required=False,default=10,
                        help="number of refinements to be used")
    parser.add_argument("--no_changing_meta",action="store_true",required=False,default=False,
                        help="do not add any meta data")

    
    args = vars(parser.parse_args())
    main(**args)
 
