import xarray as xr
import numpy as np
import pandas as pd
import enum
import subprocess
from multiprocessing import Process
import time
from netCDF4 import Dataset
from scipy.interpolate import NearestNDInterpolator
import os

def create_nc(rid,uid,lat,lon,filename):
    nc = Dataset('/hard/dongwz/CoLM-U/srftest/'+str(filename),'w',format='NETCDF4')

    # dimension
    nc.createDimension('lat',1200)
    nc.createDimension('lon',1200)

    # coordinate vars
    lats = nc.createVariable('lat',np.float64,('lat',))
    lons = nc.createVariable('lon',np.float64,('lon',))

    # var
    r_id = nc.createVariable('REGION_ID',np.int,('lat','lon'))
    u_id = nc.createVariable('URBAN_DENSITY_CLASS',np.int,('lat','lon'))

    # var atts
    lats.units = 'degrees_north'
    lons.units = 'degrees_east'
    
    lats.long_name = 'Latitude'
    lons.long_name = 'Longitude'
    
    r_id.units = 'unitless'
    u_id.units = '1-TBD, 2-HD, 3-MD'

    r_id.long_name = 'Region ID'
    u_id.long_name = 'Urban Density Class'

    # put var
    lats[:] = lat
    lons[:] = lon

    r_id[:,:] = rid
    u_id[:,:] = uid

def local_process_unit(ncarfile, modfile,cnt):
    
    nc1 = xr.open_dataset('./NCAR/'+ncarfile)
    nc2 = xr.open_dataset('./MOD/'+modfile)
    
    lat = np.array(nc1['lat'])
    lon = np.array(nc1['lon'])
    rid = np.array(nc1['REGION_ID'])
    uid = np.array(nc1['URBAN_DENSITY_CLASS'])
    mod = np.array(nc2['PCT_URBAN'])
    lc  = np.array(nc2['LC'])

    #print(rid)
    if ((~np.any(rid)) & np.any(mod)):
       print(np.any(rid))
       print('error')
       cnt=cnt+1
       return cnt

    if np.any(rid):
       i=np.where(~(rid==0))
       interp = NearestNDInterpolator(np.transpose(i), rid[i])
       filled_data = interp(*np.indices(rid.shape))
       uid = np.where(lc==0,0,uid)
       filled_data = np.where(lc==0,0,filled_data)
       create_nc(filled_data,uid,lat,lon,ncarfile)
#ncarfile = subprocess.getoutput('ls *NCAR.nc').split()#.sort()
#modfile  = subprocess.getoutput('ls *MOD2000.nc').split()#.sort()
#count = 0
ncarfile = os.listdir('./NCAR/')
#modfile  = os.listdir('./MOD/')
#print(ncarfile)
cnt = 0
for fil in ncarfile:
    #p = Process(target=local_process_unit,args=(fil,modfile[filid]))
    #p.start()
    #count=count+1
    l = len(fil)
    #print(fil[0:l-8])
    #print(os.path.splitext(fil)[0][0:l-4])
    modfile = fil[0:l-8]+'.MOD2000.nc'
    #print(modfile)
    local_process_unit(fil, modfile,cnt)
    #if count>49:
    #    p.join()
    #    count=0

print(cnt)
