
"""
Sample code for GRL submission: Capturing the Global Variability of Marine Particulate Organic Carbon Flux: A Hierarchical Bayesian Approach
Formats SST Chl-a and depth climatologies to right scales for match-up and analysis.

Created by: Pippa Edwards
"""
#%%
#packages
import xarray as xr
import numpy as np
import datetime as dt
import glob

fp = ""

#%%
#SST

#MATCH THE LONGITUDE TO CHLA AND DEPTH
#open original sst dataset
sst = xr.open_dataset(f"{fp}input_data/sst/MODEL.SST.HAD187001-198110.OI198111-202206.nc")["SST"]
#Convert lon to -180 to 180
sst = sst.assign_coords(lon=(((sst.lon + 180) % 360) - 180))
#Sort data by new longitude order
sst = sst.sortby('lon')

#remove values that are ice
sst = sst.where(sst > -1.79, np.nan)
print(np.nanmin(sst))
#save
sst.to_netcdf(f"{fp}input_data/correct_lon_SST.nc")

#make the monthly climatology
#limit to matched dataset time frame
sst = sst.sel(time = slice(dt.datetime(1997,9,1), dt.datetime(2022,7,1)))
for m in range(1,13):
    #select data that is just that month
    msst = sst.sel(time = sst.time.dt.month == m)
    msst = msst.where(msst > -1.79, np.nan)
    #make mean of those months
    msst = msst.mean(dim = "time", skipna = True)
    #assign a coordinate of month
    msst = msst.assign_coords({"month": m})
    #save
    msst.to_netcdf(f"{fp}input_data/sst/SST_{m}_monthly_climatology.nc")

#make the overall climatology
#take mean of all available data
csst = sst.mean(dim = "time", skipna= True)
csst.to_netcdf(f"{fp}input_data/sst/SST_overall_climatology.nc")

#%%
#Chla
#needs regridding to match SST 

#set slice function
#as round is confusing, use int.
#int always rounds down.
#this sets bounds for the chla data to be put into
def get_slice(l):
    l = int(l)
    if l >= 1:
        l = int(l) +1
        return l-0.0001, l-1
    elif l == 0:
        return 0, 1
    elif l < 0:
        l = int(l)
        return l, l-0.9999

#set up lat and lon
lats = np.linspace(-89.5, 89.5, 180)
lons =np.linspace(-179.5, 179.5, 360)

#do this for each year individually as files are so large.
for year in range(1997, 2023): 
    print(year)

    chla = xr.open_dataset(f"{fp}input_data/occci/original/occci_chla_{year}.nc")["chlor_a"]
    all_points = []
    for lat in lats:
        #print(lat)
        latu, latl = get_slice(lat)
        chlal = chla.sel(lat = slice(latu, latl))
        if chlal.shape[1] == 0:
            chlal = chla.sel(lat = slice(latl, latu))
        row_points = []
        for lon in lons:
            lonu, lonl = get_slice(lon)
            chlall = chlal.sel(lon = slice(lonu, lonl))
            if chlall.shape[2] == 0:
                chlall = chlal.sel(lon = slice(lonl, lonu))
            pointll = chlall.mean(dim = ["lat", "lon"], skipna = True)
            pointll = pointll.expand_dims({"lat": [lat], "lon": [lon]})
            row_points.append(pointll)
        row = xr.concat(row_points, dim="lon")
        all_points.append(row)
    new_chla = xr.concat(all_points, dim = "lat")
    new_chla.to_netcdf(f"{fp}input_data/occci/regrid/occci_chla_{year}.nc")

#%%
#create monthly climatologies of chla using same process as SST
#add all the new files into a lit

chla_files = glob.glob(f"{fp}input_data/occci/regrid/*.nc")

years = []
#make them into one big file now they are smaller (may take some time)
for f in chla_files:
    
    chlat = xr.open_dataset(f)["chlor_a"]
    years.append(chlat)
chla_all = xr.concat(years, dim = "time")
#limit to same range as SST
chla_all = chla_all.sel(time = slice(dt.datetime(1997,9,4), dt.datetime(2022,7,1)))

#use same method as SST for monthly:
for m in range(1,13):
    chlam = chla_all.sel(time = chla_all.time.dt.month == m)
    chlam = chlam.mean(dim = "time", skipna = True)
    chlam.assign_coords({"month": m})
    chlam.to_netcdf(f"{fp}input_data/occci/{m}_occci_monthly_regrid.nc")

#and overall:
chla_all = chla_all.mean(dim = "time", skipna = True)
chla_all.to_netcdf(f"{fp}input_data/occci/occci_overall_climatology.nc")

#%%
#Depth
#this is obviously not a climatology but it is still needed!
#First bathymetric map needs regridding

#open original file
bathy = xr.open_dataset(f"{fp}input_data/bathy/GLO-MFC_001_024_mask_bathy.nc")["deptho"]

#make lat and lon grid
lats = np.linspace(-89.5, 89.5, 180)
lons =np.linspace(-179.5, 179.5, 360)

#this is slow but it gives me time to write the rest while this runs!
#make a list to convert into xarray of lats and lons
all_points = []
print("generating regridded")
for lat in lats: #for each lat
    #print(lat)
    latu, latl = get_slice(lat) #find 1o area
    bathyl = bathy.sel(latitude = slice(latu, latl)) #select values within that
    #make sure it has a shape
    if bathyl.shape[0] == 0:
        #if it doesn't repeat with switched latu and latl
        bathyl = bathy.sel(latitude = slice(latl, latu))
    
    #make a list to add lon means to 
    row_points = []
    #for each set lon
    for lon in lons:
        lonu, lonl = get_slice(lon)#find 1o area
        bathyll = bathyl.sel(longitude = slice(lonu, lonl))#select values within that
         #make sure it has a shape
        if bathyll.shape[1] == 0:
            #if it doesn't repeat with switched lonu and lonl
            bathyll = bathyl.sel(longitude = slice(lonl, lonu))
        #get the mean seafloor depth using dimensions lat and lon
        pointll = bathyll.mean(dim = ["latitude", "longitude"], skipna = True)
        #give the point new dimensions of lat and lon
        pointll = pointll.expand_dims({"lat": [lat], "lon": [lon]})
        #add the point to the list of lons at that latitude
        row_points.append(pointll)
    #make this list into an xarray when all lons have been added
    row = xr.concat(row_points, dim="lon")
    #add this xarray into a list for the arrays at each lat
    all_points.append(row)
#make this into an xarray when all lats have been added
new_bathy = xr.concat(all_points, dim = "lat")
#save new seafloor map as netcdff
new_bathy.to_netcdf(f"{fp}input_data/bathy/regridded_seafloor.nc")

#make a file for each depth horizon where at all depths => x the array will be x for 100-1500m
for x in np.arange(100,1600,100):
    #print what depth is mapping
    print(x) 
    #set up blank array 
    depthmap = np.ones((180,360))
    for lat in range(len(lats)):
        for lon in range(len(lons)):
            #find the depth in the bathymetric map
            depth  = new_bathy[lat, lon]
            #if it is greater or equal to 100 then that location can be used
            if depth >= x:
                depthmap[lat,lon] = x #set depth to x
            else:
            #if not set it to nan 
                depthmap[lat,lon] = np.nan
    #make this into an xarray so it can be saved as a netcdf
    newd = xr.Dataset({f"depth":(("lat", "lon"),
                                depthmap)},
                        coords = {"lat": np.linspace(-89.5, 89.5, 180),
                                "lon": np.linspace(-179.5, 179.5, 360)})
    newd.to_netcdf(f"{fp}input_data/bathy/depth{x}_map.nc")
#%%

depthmap = np.ones((180,360))
for lat in range(len(lats)):
    for lon in range(len(lons)):
        #find the depth in the bathymetric map
        depth = new_bathy[lat, lon]
        #if it is greater or equal to 100 then that location can be used
        if depth >= 100:
            depthmap[lat,lon] = depth #set depth to x
        else:
        #if not set it to nan 
            depthmap[lat,lon] = np.nan
#make this into an xarray so it can be saved as a netcdf
newd = xr.Dataset({f"depth":(("lat", "lon"),
                            depthmap)},
                    coords = {"lat": np.linspace(-89.5, 89.5, 180),
                            "lon": np.linspace(-179.5, 179.5, 360)})
newd.to_netcdf(f"{fp}input_data/bathy/depthSF_map.nc")

depthmap = np.ones((180,360))
for lat in range(len(lats)):
    for lon in range(len(lons)):
        #find the depth in the bathymetric map
        depth = new_bathy[lat, lon]
        d2 = depth + 100
        #if it is greater or equal to 100 then that location can be used
        if d2 >= 100:
            depthmap[lat,lon] = d2 #set depth to x
        else:
        #if not set it to nan 
            depthmap[lat,lon] = np.nan
#make this into an xarray so it can be saved as a netcdf
newd = xr.Dataset({f"depth":(("lat", "lon"),
                            depthmap)},
                    coords = {"lat": np.linspace(-89.5, 89.5, 180),
                            "lon": np.linspace(-179.5, 179.5, 360)})
newd.to_netcdf(f"{fp}input_data/bathy/depthSF100_map.nc")
# %%
