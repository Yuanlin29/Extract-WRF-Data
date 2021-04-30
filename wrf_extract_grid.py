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
   # indw=np.where(wd>=360)
   # wd[indw]=wd[indw]-360
    wd=(wd%360.0)
    return(ws,wd)
# find the nearest grid
def nearest_grid(sta_lat,sta_lon,xlat,xlon):
    difflat=sta_lat-xlat
    difflon=sta_lon-xlon
    rad=np.multiply(difflat,difflat)+np.multiply(difflon,difflon)
    aa=np.where(rad==np.min(rad))
    ind=np.squeeze(np.array(aa))
    return tuple(ind)

wrf_dir='/home/mvi/WRF_script/WRF4/CERC-EPA/2018/WRF_run'
sites_dir='/home/yuawan/data/sites_info'
data_dir='/air_models/home/yuawan/data/sites/wrf'

# read information of sites#
sites_info=pd.read_csv(os.path.join(sites_dir,'MetEireann_Station_summary.csv'),sep=',',skiprows=1)
site_name=sites_info.pop('name')
site_num=sites_info.pop('Station Number')
sta_lat=sites_info.pop('Latitude')
sta_lon=sites_info.pop('Longitude')
print(np.shape(sta_lat))

# running time
start_time='2018-01-01'
end_time='2018-12-31'
freq_time='1d'
dates=pd.date_range(start=start_time,end=end_time,freq=freq_time)

method=1
var_num=7
domain1=4
domain2=4

for da in range(domain1,domain2+1):
    #read info of long & lat
    sample_name='wrfout_d0'+str(da)+'_{}-{:02d}-{:02d}_00:00:00'.format(dates[0].year,dates[0].month,dates[0].day)
    sample_path=os.path.join(wrf_dir,sample_name)
    sample_data=nc.Dataset(sample_path,'r')
    tmp_lat=sample_data.variables['XLAT']
    tmp_lon=sample_data.variables['XLONG']
    xlat=tmp_lat[0]
    xlon=tmp_lon[0]
    
    #output path
    data_path=os.path.join(data_dir,'d0'+str(da),'test1')
    os.mkdir(data_path)
    print(data_path)
    
    for isite in range(len(site_name)):
        var=np.zeros((len(dates)*24,var_num))
        time=np.zeros((len(dates)*24,4))
        iline=0
        if(method==1):
           ind=nearest_grid(sta_lon[isite],sta_lat[isite],xlon,xlat)
           print(ind)

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
                   temp=T[itime,ind[0],ind[1]]-273.15
                   psfc=P[itime,ind[0],ind[1]]
                   u10 =U10[itime,ind[0],ind[1]]
                   v10 =V10[itime,ind[0],ind[1]]
                   rainc=RAINC[itime,ind[0],ind[1]]   
                   rainnc=RAINNC[itime,ind[0],ind[1]]   
                   (ws,wd)=uv2wsd(u10,v10)
                   rain=rainc+rainnc
               
                   # print(iline) 
                   var[iline,:]=[temp,psfc,u10,v10,ws,wd,rain]
                   time[iline,:]=[date.year,date.month,date.day,itime]
                   iline=iline+1
               final=np.concatenate([time,var],axis=1)
               final_df=pd.DataFrame(final,columns=['Year','Month','Day','Hour','Temp','Pressure','U10','V10','Ws','Wd','Rain'])
               # control the format
               final_df['Year']=final_df['Year'].map(lambda x : '%.0f' %x)
               final_df['Month']=final_df['Month'].map(lambda x : '%.0f' %x)
               final_df['Day']=final_df['Day'].map(lambda x : '%.0f' %x)
               final_df['Hour']=final_df['Hour'].map(lambda x : '%.0f' %x)
           #output file for each site
           output_name=site_name[isite]+'.csv'
           print(output_name)
           output_path=os.path.join(data_path,output_name)
           final_df.to_csv(output_path,sep=',',index=False,float_format='%7.3f')   
     
