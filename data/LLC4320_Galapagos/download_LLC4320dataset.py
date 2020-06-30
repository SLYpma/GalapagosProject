import os, sys
import numpy as np
import xarray as xr
from matplotlib import pyplot as plt
import xmitgcm as xm
from netCDF4 import Dataset
from xmitgcm import llcreader
import fsspec
from datetime import datetime
from datetime import timedelta
import pandas as pd
from os import path
from glob import glob

number_of_run = int(sys.argv[1])
crash_day = int(sys.argv[2])
number_of_days = 10

datelist = pd.date_range(datetime(2011,9,13,0), periods=416).tolist()
days_to_save = np.arange(number_of_run*number_of_days+crash_day,number_of_run*number_of_days+number_of_days,1) 
data_out = "/projects/0/topios/hydrodynamic_data/LLC4320_Galapagos/hourly_output/"
time_reference = datetime(2011,9,10,0)
delta_t = 25
model = llcreader.ECCOPortalLLC4320Model()

for day in days_to_save:
    DayStr = datelist[day].strftime("%d")
    MonthStr = datelist[day].strftime("%m")
    YearStr = datelist[day].strftime("%Y")
    fname = data_out + "LLC4320_Galapagos_" + YearStr + MonthStr + DayStr + ".nc"

    start_index = (datelist[day]-time_reference).total_seconds()/delta_t
    end_index = (datelist[day+1]-time_reference).total_seconds()/delta_t
    
    ds = model.get_dataset(varnames=['U','V'], iter_start=start_index, iter_stop=end_index, k_levels=[0], read_grid = True)

    UU = np.zeros((24,1500,1500),dtype="float32")
    VV = np.zeros((24,1500,1500),dtype="float32")
    times = np.zeros((24))

    for hour in range(24):       
        time_value = start_index*delta_t + hour*3600
        # note that U=V and V=U
        V_galapagos = ds.U.isel(time=hour, k=0, face=11, i_g=slice(0,1500,1), j=slice(1000,2500,1))
        U_galapagos = ds.V.isel(time=hour, k=0, face=11, i=slice(0,1500,1), j_g=slice(1000,2500,1))

        # Get the data local
        U = np.array(U_galapagos[:,:],dtype="float32")
        V = np.array(V_galapagos[:,:],dtype="float32")
        print('loaded hour ' + str(hour) + ' on day ' + str(day))
        
        # make a new dataset to save to Netcdf
        UU[hour,:,:] = np.expand_dims(np.rot90(U),axis=0)
        VV[hour,:,:] = np.expand_dims(np.rot90(V),axis=0)
        times[hour] = np.expand_dims(time_value,axis=0)

    dsU = xr.DataArray(UU,dims=['time','YG','XC'],name="UVEL")
    dsU.attrs['standard_name']='Uvel'
    dsU.attrs['long_name']='Zonal Component of Velocity (m/s)'
    dsU.attrs['units']='m/s'
    dsU.attrs['mate']='VVEL'

    dsV = xr.DataArray(VV,dims=['time','YC','XG'],name="VVEL")
    dsV.attrs['standard_name']='Vvel'
    dsV.attrs['long_name']='Meridional Component of Velocity (m/s)'
    dsV.attrs['units']='m/s'
    dsV.attrs['mate']='UVEL'

    dsT = xr.DataArray(times,dims='time',name='time')
    dsT.attrs['units']='seconds since 2011-09-10T00:00:00'
    dsT.attrs['time_origin']='2011-09-10 00:00:00'
    dsT.attrs['calendar']='gregorian'
    dsT.attrs['axis']='T'

    ds_fields = xr.merge([dsU,dsV])
    ds_fields.coords['time']=dsT 

    #save as netcdf
    ds_fields.to_netcdf(path=fname)
    print('finished saving' + fname)
