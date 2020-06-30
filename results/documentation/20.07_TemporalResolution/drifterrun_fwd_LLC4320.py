from datetime import timedelta as delta
from os import path
from glob import glob
import numpy as np
import dask
import math
import xarray as xr
from netCDF4 import Dataset
import warnings
import matplotlib.pyplot as plt
warnings.simplefilter('ignore', category=xr.SerializationWarning)

from parcels import AdvectionRK4
from parcels import Field
from parcels import FieldSet
from parcels import JITParticle
from parcels import ParticleFile
from parcels import ParticleSet
from parcels import Variable

#input
wstokes = False          #False || True
data_in_waves = "/projects/0/topios/hydrodynamic_data"
data_in_mit = "/projects/0/topios/hydrodynamic_data/LLC4320_Galapagos"
data_out = "/home/sypmauu/GalapagosProject/results/data_output"
filename_out = "DrifterRelease_200701"
galapagos_domain = [-94, -87, -3.5, 3]
seeding_distance = 1 #unit: lon/lat degree
seeding_resolution = 4 #unit: gridpoints
seeding_frequency = 1 #unit: days
advection_duration = 10 #unit: days
output_frequency = 1 #unit: hours
length_simulation = 10 #unit: days

#Get indices for Galapagos domain to run simulation
def getclosest_ij(lats,lons,latpt,lonpt):    
    """Function to find the index of the closest point to a certain lon/lat value."""
    dist_lat = (lats-latpt)**2                      # find squared distance of every point on grid
    dist_lon = (lons-lonpt)**2
    minindex_lat = dist_lat.argmin()                # 1D index of minimum dist_sq element
    minindex_lon = dist_lon.argmin()
    return minindex_lat, minindex_lon                # Get 2D index for latvals and lonvals arrays from 1D index

dfile = Dataset(data_in_mit+'/domain/LLC4320_galapagosgrid.nc')
lon = dfile.variables['longitude'][:]
lat = dfile.variables['latitude'][:]
iy_min, ix_min = getclosest_ij(lat, lon, galapagos_domain[2], galapagos_domain[0])
iy_max, ix_max = getclosest_ij(lat, lon, galapagos_domain[3], galapagos_domain[1])

#MITgcm field
varfiles = sorted(glob(data_in_mit + "/hourly_output/LLC4320_Galapagos_*.nc"))
meshfile = glob(data_in_mit+"/domain/LLC4320_galapagosgrid.nc")
files_MITgcm = {'U': {'lon': meshfile, 'lat': meshfile, 'data': varfiles},
                'V': {'lon': meshfile, 'lat': meshfile, 'data': varfiles}}
variables_MITgcm = {'U': 'UVEL', 'V': 'VVEL'}
dimensions_MITgcm = {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}
indices_MITgcm = {'lon': range(ix_min,ix_max), 'lat': range(iy_min,iy_max)}

fieldset_MITgcm = FieldSet.from_c_grid_dataset(files_MITgcm,
                                               variables_MITgcm, 
                                               dimensions_MITgcm, 
                                               indices = indices_MITgcm,
                                               tracer_interp_method='cgrid_velocity')

#Stokes Field
if wstokes:
    files_stokes = sorted(glob(data_in_waves + "/WaveWatch3data/CFSR/WW3-GLOB-30M_20[08-12]*_uss.nc"))

    variables_stokes = {'U': 'uuss',
                        'V': 'vuss'}
    dimensions_stokes = {'U': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'},
                         'V': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}}
    indices_stokes = {'lon': range(120, 220), 'lat': range(142, 170)}
    
    fieldset_stokes = FieldSet.from_netcdf(files_stokes, 
                                           variables_stokes, 
                                           dimensions_stokes,
                                           indices=indices_stokes)
    fieldset_stokes.add_periodic_halo(zonal=True, meridional=False, halosize=5)
    fieldset = FieldSet(U=fieldset_MITgcm.U + fieldset_stokes.U,
                        V=fieldset_MITgcm.V + fieldset_stokes.V)

    fname = path.join(data_out, filename_out + "_wstokes.nc")
else:
    fieldset = fieldset_MITgcm
    fname = path.join(data_out, filename_out + ".nc")   

fU=fieldset_MITgcm.U

# get all lon, lat that are land
fieldset_MITgcm.computeTimeChunk(fU.grid.time[0], 1)
lon = np.array(fU.lon[:]) 
lat = np.array(fU.lat[:])
LandMask = fU.data[0,:,:]
LandMask = np.array(LandMask)
land = np.where(LandMask == 0)

# seed particles at seeding_distance from land
lons = np.array(fU.lon[::seeding_resolution])
lats = np.array(fU.lat[::seeding_resolution])
yy, xx = np.meshgrid(lats,lons)
xcoord = np.reshape(xx,len(lons)*len(lats))
ycoord = np.reshape(yy,len(lons)*len(lats))

startlon=[]
startlat=[]

for i in range(xcoord.shape[0]):
    dist = (xcoord[i]-lon[land[1]])**2 + (ycoord[i]-lat[land[0]])**2
    minindex = dist.argmin()
    if dist[minindex]<seeding_distance and dist[minindex] != 0:
        startlon.append(xcoord[i])
        startlat.append(ycoord[i])
        
#functions to add to the kernel
def Age(fieldset, particle, time):
    particle.age = particle.age + math.fabs(particle.dt)
    if particle.age > 10*86400:
        particle.delete()
    
#additional features of the particles        
class GalapagosParticle(JITParticle):
    age = Variable('age', initial = 0.)   

def DeleteParticle(particle, fieldset, time):
    particle.delete()
    
# set particle conditions
pset = ParticleSet(fieldset=fieldset,
                   pclass=GalapagosParticle,
                   lon=startlon,
                   lat=startlat,
                   repeatdt=delta(days=seeding_frequency))

outfile = pset.ParticleFile(name=fname, outputdt=delta(hours=output_frequency))
kernels = pset.Kernel(AdvectionRK4) + pset.Kernel(Age)

pset.execute(kernels,
             runtime=delta(days=length_simulation),
             dt=delta(minutes=20),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

pset.repeatdt = None

pset.execute(kernels,
             runtime=delta(days=advection_duration),
             dt=delta(minutes=20),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

outfile.export()
outfile.close()  
        
