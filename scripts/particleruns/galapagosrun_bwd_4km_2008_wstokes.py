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
chunk_mode = False #False || 'auto' || 'specific'
wstokes = True          #False || True

if chunk_mode == 'auto':
    dask.config.set({'array.chunk-size': '2MiB'})
else:
    dask.config.set({'array.chunk-size': '128MiB'})

    
#MITgcm field
data_path_mitgcm = path.join("/home/sypmauu/GalapagosProject/data", "MIT4km/")
data_file_mitgcm = path.join(data_path_mitgcm, "RGEMS3_2008_Surf.nc")
mesh_mask_mitgcm = path.join(data_path_mitgcm, "RGEMS3_Surf_grid.nc")

filenames_mitgcm = {'U': {'lon': mesh_mask_mitgcm, 'lat': mesh_mask_mitgcm, 'data': data_file_mitgcm},
                    'V': {'lon': mesh_mask_mitgcm, 'lat': mesh_mask_mitgcm, 'data': data_file_mitgcm}}
variables_mitgcm = {'U': 'UVEL',
                    'V': 'VVEL'}
dimensions_mitgcm = {'U': {'lon': 'XG', 'lat': 'YG', 'time': 'time'},
                     'V': {'lon': 'XG', 'lat': 'YG', 'time': 'time'}}
chs_mitgcm = False
if chunk_mode == 'auto':
    chs_mitgcm = 'auto'
elif chunk_mode == 'specific':
    chs_mitgcm = {'U': {'time': 1, 'XC': 80, 'XG': 80, 'YC': 60, 'YG': 60},
                  'V': {'time': 1, 'XC': 80, 'XG': 80, 'YC': 60, 'YG': 60}}
fieldset_MITgcm = FieldSet.from_c_grid_dataset(filenames_mitgcm, variables_mitgcm, dimensions_mitgcm, 
                                               field_chunksize=chs_mitgcm, 
                                               time_periodic=delta(days=366), 
                                               tracer_interp_method='cgrid_velocity')


#Stokes Field
if wstokes:
    data_path_stokes = path.join("/projects/0/topios/hydrodynamic_data", "WaveWatch3data", "CFSR")
    data_files_stokes = sorted(glob(data_path_stokes + "/WW3-GLOB-30M_2008*_uss.nc"))
    
    variables_stokes = {'U': 'uuss',
                        'V': 'vuss'}
    dimensions_stokes = {'U': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'},
                         'V': {'lon': 'longitude', 'lat': 'latitude', 'time': 'time'}}
    chs_stokes = False
    if chunk_mode == 'auto':
        chs_stokes = 'auto'
    elif chunk_mode == 'specific':
        chs_stokes = {'time': 1, 'latitude': 32, 'longitude': 16}
    
    fieldset_stokes = FieldSet.from_netcdf(data_files_stokes, variables_stokes, dimensions_stokes, 
                                           field_chunksize=chs_stokes, 
                                           time_periodic=delta(days=366))
    fieldset_stokes.add_periodic_halo(zonal=True, meridional=False, halosize=5)
    fieldset = FieldSet(U=fieldset_MITgcm.U + fieldset_stokes.U, 
                        V=fieldset_MITgcm.V + fieldset_stokes.V)
    
#    fU = fieldset.U[0]
    fname = "/home/sypmauu/GalapagosProject/results/data_output/galapagosparticles_bwd_4km_2008_wstokes.nc"
else:
    fieldset = fieldset_MITgcm
#    fU = fieldset.U
    fname = "/home/sypmauu/GalapagosProject/results/data_output/galapagosparticles_bwd_4km_2008.nc"

#initialize where to start particles
galapagos_extent = [-91.8, -89, -1.4, 0.7]
startlon, startlat = np.meshgrid(np.arange(galapagos_extent[0],
                                           galapagos_extent[1],0.2),
                                 np.arange(galapagos_extent[2],
                                           galapagos_extent[3],0.2))

#functions to add to the kernel
def Age(fieldset, particle, time):
    particle.age = particle.age + math.fabs(particle.dt)
    if particle.age > 300*86400:
        particle.delete()      

def DeleteParticle(particle, fieldset, time):
    particle.delete()            
    
#additional features of the particles        
class GalapagosParticle(JITParticle):
    age = Variable('age', initial = 0.)    
    
# set particle conditions
pset = ParticleSet(fieldset=fieldset,
                   pclass=GalapagosParticle,
                   lon=startlon,
                   lat=startlat)
#                   time=fU.grid.time[-1])
#                   repeatdt=delta(days=5))

outfile = pset.ParticleFile(name=fname, outputdt=delta(days=1))
kernels = pset.Kernel(AdvectionRK4) + pset.Kernel(Age)  

pset.execute(kernels,
             runtime=delta(days=30),
             dt=delta(hours=-1),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

#pset.repeatdt = None

#pset.execute(AdvectionRK4+pset.Kernel(Age),
#             runtime=delta(days=300),
#             dt=delta(hours=-1),
#             output_file=outfile,
#             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

outfile.export()
outfile.close()
