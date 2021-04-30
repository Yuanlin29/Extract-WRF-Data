import numpy as np
import netCDF4 as nc
import wrf
import pandas as pd
import os
from scipy.interpolate import griddata

# convert uv to ws,wd
def uv2wsd(x,y):
    ws=np.sqrt(x**2+y**2)
    wd=np.arctan2(y,x)
    wd=wd*180/np.pi
    wd=270-wd
    ind=np.where(wd>=360)
    wd[ind]=wd[ind]-360
    return(ws,wd)

wrf_dir='/home/mvi/WRF_script/WRF4/CERC-EPA/2018/WRF_run'
sites_dir='/home/yuawan/data/sites_info'
data_dir='/air_models/home/yuawan/data/sites/wrf'

# read information of sites#
sites_info=pd.read_csv(os.path.join(sites_dir,'MetEireann_Station_summary.csv'),sep=',',skiprows=1)
site_name=sites_info.pop('name')
lat=sites_info.pop('Latitude')
lon=sites_info.pop('Longitude')

# running time
start_time='2018-01-01'
end_time='2018-12-31'
freq_time='1d'
dates=pd.date_range(start=start_time,end=end_time,freq=freq_time)

var_num=7
domain1=4
domain2=4

for da in range(domain1,domain2+1):
    #output path
    data_path=os.path.join(data_dir,'d0'+str(da),'test2')
    os.mkdir(data_path)
    print(data_path)

    var=np.zeros((len(dates)*24,len(site_name),var_num))
    time=np.zeros((len(dates)*24,4))
    iline=0

    # read LON & LAT from wrfout
    sample_name='wrfout_d0'+str(da)+'_{}-{:02d}-{:02d}_00:00:00'.format(dates[0].year,dates[0].month,dates[0].day)
    sample_path=os.path.join(wrf_dir,sample_name)
    sample_data=nc.Dataset(sample_path,'r')
    tmp_lat=sample_data.variables['XLAT']
    tmp_lon=sample_data.variables['XLONG']
    LAT=tmp_lat[0]
    LON=tmp_lon[0]
    points=np.c_[LON.flatten().T,LAT.flatten().T]

    for date in dates:    
        print(date)

        file_name='wrfout_d0'+str(da)+'_{}-{:02d}-{:02d}_00:00:00'.format(date.year,date.month,date.day)
        file_path=os.path.join(wrf_dir,file_name)
        wrfin=nc.Dataset(file_path,'r')
        
        # variables
        T=wrfin.variables['T2']
        P=wrfin.variables['PSFC']
        U10=wrfin.variables['U10']
        V10=wrfin.variables['V10']
        RAINC=wrfin.variables['RAINC']
        RAINNC=wrfin.variables['RAINNC']
    
        for itime in range(0,24):
            temp=griddata(points,T[itime,:,:].flatten(),(lon,lat),method='linear')
            psfc=griddata(points,P[itime,:,:].flatten(),(lon,lat),method='linear')
            u10 =griddata(points,U10[itime,:,:].flatten(),(lon,lat),method='linear') 
            v10 =griddata(points,V10[itime,:,:].flatten(),(lon,lat),method='linear')
            rainc=griddata(points,RAINC[itime,:,:].flatten(),(lon,lat),method='linear')
            rainnc=griddata(points,RAINNC[itime,:,:].flatten(),(lon,lat),method='linear')
            
            temp=temp-273.15
            (ws,wd)=uv2wsd(u10,v10)
            rain=rainc+rainnc

            tmp=(np.c_[temp.flatten().T,psfc.flatten().T,u10.flatten().T,v10.flatten().T,ws.flatten().T,wd.flatten().T,rain.flatten().T])
            #different arraies
            var[iline,:,:]=tmp
            time[iline,:]=[date.year,date.month,date.day,itime]
            iline=iline+1
   
    # output file for each site
    for isite in range(len(site_name)):
        final=np.concatenate([time,var[:,isite,:]],axis=1)
        final_df=pd.DataFrame(final,columns=['Year','Month','Day','Hour','Temp','Pressure','U10','V10','Ws','Wd','Rain'])
        # control the format
        final_df['Year']=final_df['Year'].map(lambda x : '%.0f' %x)
        final_df['Month']=final_df['Month'].map(lambda x : '%.0f' %x)
        final_df['Day']=final_df['Day'].map(lambda x : '%.0f' %x)
        final_df['Hour']=final_df['Hour'].map(lambda x : '%.0f' %x)
        output_name=site_name[isite]+'.csv'
        print(output_name)
        output_path=os.path.join(data_path,output_name)
        final_df.to_csv(output_path,sep=',',index=False,float_format='%7.3f')


