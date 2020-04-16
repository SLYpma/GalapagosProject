from parcels import FieldSet, Field, ParticleSet, JITParticle, AdvectionRK4, ErrorCode, Variable
from datetime import timedelta as delta
from glob import glob
import numpy as np
import xarray as xr
import os
import warnings
warnings.simplefilter('ignore', category=xr.SerializationWarning)

ddir = "/home/sypmauu/GalapagosProject/data/MIT4km/"
wstokes = True

# set field
varfiles = glob(ddir+'RGEMS3_2008_Surf.nc')
meshfile = glob(ddir+'RGEMS3_Surf_grid.nc')
MITgcm_files = {'U': {'lon': meshfile, 'lat': meshfile, 'data': varfiles},
                'V': {'lon': meshfile, 'lat': meshfile, 'data': varfiles}}
MITgcm_variables = {'U': 'UVEL', 'V': 'VVEL'}
MITgcm_dimensions = {'lon': 'XG', 'lat': 'YG', 'time': 'time'}
fieldset_MITgcm = FieldSet.from_c_grid_dataset(MITgcm_files, MITgcm_variables, 
                                               MITgcm_dimensions,
#                                               time_periodic=delta(days=366), 
                                               tracer_interp_method='cgrid_velocity')

if wstokes:
    stokesfiles = sorted(glob('/projects/0/topios/hydrodynamic_data/WaveWatch3data/CFSR/WW3-GLOB-30M_2008*_uss.nc'))
    stokesdimensions = {'lat': 'latitude', 'lon': 'longitude', 'time': 'time'}
    stokesvariables = {'U': 'uuss', 'V': 'vuss'}
    fieldset_stokes = FieldSet.from_netcdf(stokesfiles, stokesvariables, 
                                           stokesdimensions) 
#                                           time_periodic=delta(days=366))
    fieldset_stokes.add_periodic_halo(zonal=True, meridional=False, halosize=5)

    fieldset = FieldSet(U=fieldset_MITgcm.U+fieldset_stokes.U, 
                        V=fieldset_MITgcm.V+fieldset_stokes.V)
    fU = fieldset.U[0]
    fname = "/home/sypmauu/GalapagosProject/results/data_output/galapagosparticles_bwd_4km_2008_wstokes.nc"
else:
    fieldset = fieldset_MITgcm
    fU = fieldset.U
    fname = "/home/sypmauu/GalapagosProject/results/data_output/galapagosparticles_bwd_4km_2008.nc"

fieldset.computeTimeChunk(fU.grid.time[-1], -1)

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
                   lat=startlat,
                   time=fU.grid.time[-1],
                   repeatdt=delta(days=1))

outfile = pset.ParticleFile(name=fname, outputdt=delta(days=1))

pset.execute(AdvectionRK4+pset.Kernel(Age),
             runtime=delta(days=30),
             dt=delta(hours=-1),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

pset.repeatdt = None

pset.execute(AdvectionRK4+pset.Kernel(Age),
             runtime=delta(days=270),
             dt=delta(hours=-1),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

outfile.export()
outfile.close()
