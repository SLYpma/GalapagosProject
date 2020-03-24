import cartopy
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

def PlotTrajectories(plon,plat,map_extent=[-110, -70, -20, 30]):
    
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
        
        
def ReadTrajectories(namefile):
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