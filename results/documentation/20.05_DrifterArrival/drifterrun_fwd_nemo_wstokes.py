from datetime import timedelta as delta
from os import path
from glob import glob
import numpy as np
import dask
import math
import xarray as xr
import warnings
warnings.simplefilter('ignore', category=xr.SerializationWarning)

from parcels import AdvectionRK4
from parcels import Field
from parcels import FieldSet
from parcels import JITParticle
from parcels import ParticleFile
from parcels import ParticleSet
from parcels import Variable

#input
wstokes = True          #False || True
data_in = "/projects/0/topios/hydrodynamic_data"
data_out = "/home/sypmauu/GalapagosProject/results/data_output"
filename_out = "DrifterArrival_fwd_nemo"


#NEMO field
ufiles = sorted(glob(data_in + "/NEMO-MEDUSA/ORCA0083-N006/means/ORCA0083-N06_200[8-9]*d05U.nc"))
vfiles = [u.replace('05U.nc', '05V.nc') for u in ufiles]
meshfile = glob(data_in + "/NEMO-MEDUSA/ORCA0083-N006/domain/coordinates.nc")

files_nemo = {'U': {'lon': meshfile, 'lat': meshfile, 'data': ufiles},
              'V': {'lon': meshfile, 'lat': meshfile, 'data': vfiles}}
variables_nemo = {'U': 'uo', 'V': 'vo'}
dimensions_nemo = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}
indices_nemo = {'lon': range(2005, 2605), 'lat': range(1410, 1578)}

fieldset_nemo = FieldSet.from_nemo(files_nemo, 
                                   variables_nemo,
                                   dimensions_nemo,
                                   indices=indices_nemo)


#Stokes Field
if wstokes:
    files_stokes = sorted(glob(data_in + "/WaveWatch3data/CFSR/WW3-GLOB-30M_200[8-9]*_uss.nc"))

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
    fieldset = FieldSet(U=fieldset_nemo.U + fieldset_stokes.U,
                        V=fieldset_nemo.V + fieldset_stokes.V)

    fU = fieldset.U[0]
    fname = path.join(data_out, filename_out + "_wstokes.nc")
else:
    fieldset = fieldset_nemo
    fU = fieldset.U
    fname = path.join(data_out, filename_out + ".nc")

    
#initialize where to start particles
GC = [-90.5,-0.5] #centre of the Galapagos islands
release_extent = [GC[0]-3.5, GC[0]+3.5, GC[1]-3.5, GC[1]+3.5]
width = 7 #thicknes (in number of particles) of release square around Galapagos
templon, templat = np.meshgrid(np.arange(release_extent[0],
                                         release_extent[1],0.2),
                               np.arange(release_extent[2],
                                         release_extent[3],0.2))
templon[width:-width,width:-width]= np.nan
templat[width:-width,width:-width]= np.nan
startlon = templon[np.logical_not(np.isnan(templon))]
startlat = templat[np.logical_not(np.isnan(templat))]


#mask for different islands galapagos 
lenx = len(fU.grid.lon)
leny = len(fU.grid.lat)
galapagosmask = np.zeros((leny,lenx))
extent1 = [-91.7, -89.9, -1.1, 0.2]
extent2 = [-90.8, -90.3, 0.2, 0.7]
extent3 = [-90.6, -89.5, -1.6, -1.1]
extent4 = [-89.9, -89, -1.1, -0.5]
for x in range(0, lenx):
    for y in range(0, leny):
        if (fU.grid.lon[x] >= extent1[0] and 
            fU.grid.lon[x] < extent1[1] and
            fU.grid.lat[y] >= extent1[2] and 
            fU.grid.lat[y] < extent1[3]):
            galapagosmask[y, x] = 1
        if (fU.grid.lon[x] >= extent2[0] and 
            fU.grid.lon[x] < extent2[1] and
            fU.grid.lat[y] >= extent2[2] and 
            fU.grid.lat[y] < extent2[3]):
            galapagosmask[y, x] = 2
        if (fU.grid.lon[x] >= extent3[0] and 
            fU.grid.lon[x] < extent3[1] and
            fU.grid.lat[y] >= extent3[2] and 
            fU.grid.lat[y] < extent3[3]):
            galapagosmask[y, x] = 3
        if (fU.grid.lon[x] >= extent4[0] and 
            fU.grid.lon[x] < extent4[1] and
            fU.grid.lat[y] >= extent4[2] and 
            fU.grid.lat[y] < extent4[3]):
            galapagosmask[y, x] = 4
lon = fU.grid.lon
lat = fU.grid.lat     
depth = 0
fieldset.add_field(Field('galapagosmask', galapagosmask, 
                         lon=lon,lat=lat,depth=depth,
                         mesh='spherical', interp_method='nearest',
                         allow_time_extrapolation=True))


#functions to add to the kernel
def Age(fieldset, particle, time):
    particle.age = particle.age + math.fabs(particle.dt)
    if particle.age > 180*86400:
        particle.delete()

def DeleteParticle(particle, fieldset, time):
    particle.delete()

def SampleGalapagos(fieldset, particle, time):
    if fieldset.galapagosmask[time, particle.depth, particle.lat, particle.lon] == 1:
        particle.visitedgalapagos = 1
    if fieldset.galapagosmask[time, particle.depth, particle.lat, particle.lon] == 2:
        particle.visitedgalapagos = 2
    if fieldset.galapagosmask[time, particle.depth, particle.lat, particle.lon] == 3:
        particle.visitedgalapagos = 3
    if fieldset.galapagosmask[time, particle.depth, particle.lat, particle.lon] == 4:
        particle.visitedgalapagos = 4
   
    
#additional features of the particles        
class GalapagosParticle(JITParticle):
    visitedgalapagos = Variable('visitedgalapagos', initial=0.)
    age = Variable('age', initial = 0.)   
    
    
# set particle conditions
pset = ParticleSet(fieldset=fieldset,
                   pclass=GalapagosParticle,
                   lon=startlon,
                   lat=startlat,
                   repeatdt=delta(days=5))

outfile = pset.ParticleFile(name=fname, outputdt=delta(days=1))
kernels = pset.Kernel(AdvectionRK4) + pset.Kernel(Age) + pset.Kernel(SampleGalapagos)

pset.execute(kernels,
             runtime=delta(days=366),
             dt=delta(hours=1),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

pset.repeatdt = None

pset.execute(kernels,
             runtime=delta(days=180),
             dt=delta(hours=1),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

outfile.export()
outfile.close()