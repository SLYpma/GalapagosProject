from parcels import FieldSet, Field, ParticleSet, JITParticle, AdvectionRK4, ErrorCode, Variable
from datetime import timedelta as delta
from glob import glob
import numpy as np
import xarray as xr
import os
import warnings
warnings.simplefilter('ignore', category=xr.SerializationWarning)

ddir = "/home/sypmauu/GalapagosProject/data/MIT4km/"
fname = "/home/sypmauu/GalapagosProject/results/data_output/galapagosparticles_bwd_4km_2008_website.nc"

# set field
varfiles = glob(ddir+'RGEMS3_2008_Surf.nc')
meshfile = glob(ddir+'RGEMS3_Surf_grid.nc')
MITgcm_files = {'U': {'lon': meshfile, 'lat': meshfile, 'data': varfiles},
                'V': {'lon': meshfile, 'lat': meshfile, 'data': varfiles}}
MITgcm_variables = {'U': 'UVEL', 'V': 'VVEL'}
MITgcm_dimensions = {'lon': 'XG', 'lat': 'YG', 'time': 'time'}
fieldset = FieldSet.from_c_grid_dataset(MITgcm_files, MITgcm_variables, 
                                        MITgcm_dimensions, 
                                        tracer_interp_method='cgrid_velocity')

fU=fieldset.U
fieldset.computeTimeChunk(fU.grid.time[-1], -1)

#initialize where to start particles
galapagos_extent = [-91.8, -89, -1.4, 0.7]
startlon, startlat = np.meshgrid(np.arange(galapagos_extent[0],
                                           galapagos_extent[1],0.1),
                                 np.arange(galapagos_extent[2],
                                           galapagos_extent[3],0.1))

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
#pset = ParticleSet(fieldset=fieldset,
#                   pclass=GalapagosParticle,
#                   lon=startlon,
#                   lat=startlat,
#                   time=fU.grid.time[-1],
#                   repeatdt=delta(days=1))

pset = ParticleSet(fieldset=fieldset,
                   pclass=GalapagosParticle,
                   lon=startlon,
                   lat=startlat,
                   time=fU.grid.time[-1])

outfile = pset.ParticleFile(name=fname, outputdt=delta(hours=3))

pset.execute(AdvectionRK4+pset.Kernel(Age),
             runtime=delta(days=30),
             dt=delta(hours=-1),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

#pset.repeatdt = None

#pset.execute(AdvectionRK4+pset.Kernel(Age),
#             runtime=delta(days=270),
#             dt=delta(hours=-1),
#             output_file=outfile,
#             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

outfile.export()
outfile.close()
