import cartopy
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from glob import glob
from netCDF4 import Dataset

############## PLOTTING ##############

def PlotTrajectories(plon,plat,map_extent=[-110, -70, -20, 30]):
    """
    Plots trajectories specified by plon and plat in red
    map_extent = [minlon,maxlon,minlat,maxlat] determines crop 
    """
    projection = cartopy.crs.PlateCarree(central_longitude=0)
    fig, ax = plt.subplots(subplot_kw={'projection': projection}, figsize=(7,7))
    
    grd = ax.gridlines(draw_labels=True,
                          color='gray', 
                          alpha=0.5, 
                          linestyle='--')
    grd.xlabels_top = False
    grd.ylabels_right = False
    grd.xlabel_style = {'size': 15, 'color': 'black'}  
    grd.ylabel_style = {'size': 15, 'color': 'black'}
    ax.coastlines()
    ax.add_feature(cartopy.feature.LAND, facecolor=(160/255, 160/255, 160/255))
    ax.set_extent(map_extent) 
    
    for particle in range (len(plon)):
        ax.plot(plon[particle,:], plat[particle,:], 
                   linewidth=1,
                   color='r') 
       

    
def PlotBathyTraj(plon,plat,lon,lat,bathy,
                  figsize=(10,4),
                  map_extent=[-110, -70, -20, 30],
                  Zmin=10, Zmax=4500):
    """
    Plots trajectories specified by plon and plat in red with bathymetry in background
    map_extent = [minlon,maxlon,minlat,maxlat] determines crop 
    """
    fig, ax = plt.subplots(figsize=figsize)
        
    levels = np.linspace(Zmin, Zmax, 41)
    fig = ax.contourf(lon,lat,bathy,       
                      levels = levels,    
                      cmap='bone',       
                      extend='both',      
                      origin='lower')     
    ax.set_title('bathymetry')          
    ax.set_xlabel('longitude')                 
    ax.set_ylabel('latitude') 
    ax.set_xlim(map_extent[0:2])
    ax.set_ylim(map_extent[2:4])
    ax.set_facecolor('gray')
    cbar = plt.colorbar(fig, ax=ax)
    cbar.ax.set_ylabel('depth (m)')
    
    for particle in range (len(plon)):
        ax.plot(plon[particle,:], plat[particle,:], 
                   linewidth=1,
                   color='r')     
    
    
    
def PlotSSTquiver(lon,lat,Tvel,Uvel,Vvel,
                  figsize=(10, 8),
                  map_extent=[-120, -70, -10, 10],
                  Tmin=20, Tmax=30,                                #colorbar limits
                  dq=12, headwidth=5, headlength=7, scale=0.7):    #arrow specifics
    """
    contourf plot of SST with quiver of velocity
    map_extent = [minlon,maxlon,minlat,maxlat] determines crop
    colorbar limits can be given by Tmin and Tmax
    dq is stepsize of arrows = 12 (10 per degree in case of 1/12 degree output)
    scale changes length of arrows (smaller number, larger arrow) 
    """
    fig, ax = plt.subplots(figsize=figsize)

    levels = np.linspace(Tmin, Tmax, 41)
    fig = ax.contourf(lon,lat,Tvel,       
                levels = levels,    
                cmap='Spectral_r',       
                extend='both',      
                origin='lower')     
    ax.set_title('sea surface temperature')          
    ax.set_xlabel('longitude')                 
    ax.set_ylabel('latitude') 
    ax.set_xlim(map_extent[0:2])
    ax.set_ylim(map_extent[2:4])
    ax.set_facecolor('gray')
    cbar = plt.colorbar(fig, ax=ax)
    cbar.ax.set_ylabel('SST (\xb0C)')

    if lat.ndim == 1:
        ax.quiver(lon[0::dq],lat[0::dq],
                  Uvel[0::dq,0::dq],Vvel[0::dq,0::dq],
                  headwidth=headwidth,
                  headlength=headlength,
                  scale_units='xy',
                  angles='xy',
                  scale=scale)
    
    elif lat.ndim == 2:
        ax.quiver(lon[0::dq,0::dq],lat[0::dq,0::dq],
                  Uvel[0::dq,0::dq],Vvel[0::dq,0::dq],
                  headwidth=headwidth,
                  headlength=headlength,
                  scale_units='xy',
                  angles='xy',
                  scale=scale)    
    
def PlotBathyquiver(lon,lat,bathy,Uvel,Vvel,
                    figsize=(10, 8),
                    map_extent=[-120, -70, -10, 10],
                    Zmin=10, Zmax=4500,                              #colorbar limits
                    dq=12, headwidth=5, headlength=7, scale=0.7):    #arrow specifics
    """
    contourf plot of bathymetry with quiver of velocity
    map_extent = [minlon,maxlon,minlat,maxlat] determines crop
    colorbar limits can be given by Zmin and Zmax
    dq is stepsize of arrows = 12 (10 per degree in case of 1/12 degree output)
    scale changes length of arrows (smaller number, larger arrow) 
    """
    fig, ax = plt.subplots(figsize=figsize)

    levels = np.linspace(Zmin, Zmax, 41)
    fig = ax.contourf(lon,lat,bathy,       
                levels = levels,    
                cmap='ocean',       
                extend='both',      
                origin='lower')     
    ax.set_title('bathymetry')          
    ax.set_xlabel('longitude')                 
    ax.set_ylabel('latitude') 
    ax.set_xlim(map_extent[0:2])
    ax.set_ylim(map_extent[2:4])
    ax.set_facecolor('gray')
    cbar = plt.colorbar(fig, ax=ax)
    cbar.ax.set_ylabel('depth (m)')

    if lat.ndim == 1:
        ax.quiver(lon[0::dq],lat[0::dq],
                  Uvel[0::dq,0::dq],Vvel[0::dq,0::dq],
                  headwidth=headwidth,
                  headlength=headlength,
                  scale_units='xy',
                  angles='xy',
                  scale=scale)
    
    elif lat.ndim == 2:
        ax.quiver(lon[0::dq,0::dq],lat[0::dq,0::dq],
                  Uvel[0::dq,0::dq],Vvel[0::dq,0::dq],
                  headwidth=headwidth,
                  headlength=headlength,
                  scale_units='xy',
                  angles='xy',
                  scale=scale)     
    
############## READING ##############    
        
def ReadTrajectories(namefile):
    """
    Reads any particle .nc file (output from oceanparcels)
    new diagnostics can be added using the hasattr function, see below 
    """
    Traj = {}
    pfile = xr.open_dataset(namefile, decode_cf=True)
    Traj['lon'] = np.ma.filled(pfile.variables['lon'], np.nan)
    Traj['lat'] = np.ma.filled(pfile.variables['lat'], np.nan)
    Traj['time'] = np.ma.filled(pfile.variables['time'], np.nan)
    
    if hasattr(pfile,'age'):
        Traj['age'] = np.ma.filled(pfile.variables['age'], np.nan)
    
    if hasattr(pfile,'visitedgalapagos'):
        Traj['visitedgalapagos'] = np.ma.filled(pfile.variables['visitedgalapagos'], np.nan) 
    
    return Traj


def ReadNemo(ddir, ddate, map_crop=[2000, 3001, 1000, 2001]):
    """
    Reads lat,lon,SST,Usurface,Vsurface of one Nemofile 
    Nemofile specified by ddate
    map_crop specifies which data selection to take
    Can be extended in the future to also load other diagnostics/vertical levels etc.
    """
    Field = {}
    latdim = np.arange(map_crop[2],map_crop[3])
    londim = np.arange(map_crop[0],map_crop[1])

    Ufile = ddir+'ORCA0083-N06_'+str(ddate)+'d05U.nc'
    Vfile = ddir+'ORCA0083-N06_'+str(ddate)+'d05V.nc'
    Tfile = ddir+'ORCA0083-N06_'+str(ddate)+'d05T.nc'

    dfile = Dataset(Ufile)
    Field['Uvel'] = dfile.variables['uo'][0,0,latdim,londim]
    Field['lat'] = dfile.variables['nav_lat'][latdim,londim]
    Field['lon'] = dfile.variables['nav_lon'][latdim,londim]
    dfile = Dataset(Vfile)
    Field['Vvel'] = dfile.variables['vo'][0,0,latdim,londim]
    dfile = Dataset(Tfile)
    Field['Tvel'] = dfile.variables['sst'][0,latdim,londim]
    
    return Field


def ReadNemoMean(ddir, ddates, nz, map_crop=[2000, 3001, 1000, 2001]):
    """
    Reads lat,lon,SST,Usurface,Vsurface of Nemofiles and computes mean 
    ddates specifies range of files to be read
    nz specifies at which depth layer files should be read
    map_crop specifies which data selection to take
    Can be extended in the future to also load other diagnostics/vertical levels etc.
    """
    
    Field = {}
    latdim = np.arange(map_crop[2],map_crop[3])
    londim = np.arange(map_crop[0],map_crop[1])

    Ufiles = sorted(glob(ddir+'ORCA0083-N06_20[00-10]*d05U.nc'))
    Vfiles = sorted(glob(ddir+'ORCA0083-N06_20[00-10]*d05V.nc'))
    Tfiles = sorted(glob(ddir+'ORCA0083-N06_20[00-10]*d05T.nc'))

    Umean = np.zeros([len(latdim),len(londim)], dtype=float)
    Vmean = np.zeros([len(latdim),len(londim)], dtype=float)
    Tmean = np.zeros([len(latdim),len(londim)], dtype=float)
    
    teller = 0
    for t in ddates:
        teller += 1
        dfile = Dataset(Ufiles[t])
        Uvel = dfile.variables['uo'][0,nz,latdim,londim]
        Umean = Umean + Uvel
        dfile = Dataset(Vfiles[t])
        Vvel = dfile.variables['vo'][0,nz,latdim,londim]
        Vmean = Vmean + Vvel
        dfile = Dataset(Tfiles[t])
        Tvel = dfile.variables['sst'][0,latdim,londim]
        Tmean = Tmean + Tvel
        
    Field['Tmean'] = Tmean/teller
    Field['Umean'] = Umean/teller
    Field['Vmean'] = Vmean/teller    
    Field['lat'] = dfile.variables['nav_lat'][latdim,londim]
    Field['lon'] = dfile.variables['nav_lon'][latdim,londim]
    
    return Field
    
