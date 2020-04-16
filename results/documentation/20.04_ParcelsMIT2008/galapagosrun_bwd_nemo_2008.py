from parcels import FieldSet, Field, ParticleSet, JITParticle, AdvectionRK4, ErrorCode, Variable
from datetime import timedelta as delta
from glob import glob
import numpy as np
import xarray as xr
import os
import warnings
warnings.simplefilter('ignore', category=xr.SerializationWarning)

ddir = "/projects/0/topios/hydrodynamic_data/NEMO-MEDUSA/ORCA0083-N006/"
fname = "/home/sypmauu/GalapagosProject/results/data_output/galapagosparticles_bwd_nemo_2008.nc"

# set field, only year 2008
ufiles = sorted(glob(ddir+'means/ORCA0083-N06_2008*d05U.nc'))
vfiles = [u.replace('05U.nc', '05V.nc') for u in ufiles]
meshfile = glob(ddir+'domain/coordinates.nc')
nemofiles = {'U': {'lon': meshfile, 'lat': meshfile, 'data': ufiles},
             'V': {'lon': meshfile, 'lat': meshfile, 'data': vfiles}}
nemovariables = {'U': 'uo', 'V': 'vo'}
nemodimensions = {'lon': 'glamf', 'lat': 'gphif', 'time': 'time_counter'}
indices = {'lon': range(2005, 2605), 'lat': range(1410, 1578)}
fieldset = FieldSet.from_nemo(nemofiles, nemovariables, nemodimensions,
                              time_periodic=delta(days=366),
                              indices=indices)

fU = fieldset.U
#fieldset.computeTimeChunk(fU.grid.time[0], 1)    #forward
fieldset.computeTimeChunk(fU.grid.time[-1], -1)  #backward

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
                   repeatdt=delta(days=5))

outfile = pset.ParticleFile(name=fname, outputdt=delta(days=1))

pset.execute(AdvectionRK4+pset.Kernel(Age),
             runtime=delta(days=365),
             dt=delta(hours=-1),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

pset.repeatdt = None

pset.execute(AdvectionRK4+pset.Kernel(Age),
             runtime=delta(days=300),
             dt=delta(hours=-1),
             output_file=outfile,
             recovery={ErrorCode.ErrorOutOfBounds: DeleteParticle})

outfile.export()
outfile.close()
