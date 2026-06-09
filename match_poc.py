"""
Sample code for GRL submission: Capturing the Global Variability of Marine Particulate Organic Carbon Flux: A Hierarchical Bayesian Approach
Code makes match between POC flux, OC-CCI and SST dataset in time and space
Plots figure S1

Created by: Pippa Edwards
"""

#%%
#Make the matches between the data frames

#import packages
import pandas as pd
import numpy as np
import xarray as xr
import datetime as dt
import warnings
warnings.filterwarnings("ignore")
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

fp = ""
#%%

#set up small function for setting latitude and longitudes for SST match
def temp_match(l):
    l = int(l)
    if l >= 0:
        l = l + 0.5
    else:
        l = l - 0.5
    return l

#load POC base database depending on which one is being used
alldf = pd.read_csv(f"{fp}/Global_POC_Database_2026-03-10.csv")

print(alldf.columns)
"""
Index(['Latitude [°]', 'Longitude [°]', 'Depth [m]', 'Lat_grid [°]',
       'Lon_grid [°]', 'Depth_grid [m]', 'Land_flag', 'Ocean', 'Elevation',
       'Date', 'Date_deployment', 'Date_recovery', 'Duration', 'Duration_unit',
       'Date_num', 'Year', 'Month', 'Season', 'POC_converted [mg/m-2/d-1]',
       'POC_method', 'POC_comment', 'Instrument', 'General_category',
       'ID_Publisher', 'Reference', 'Investigator', 'POC_raw', 'POC_raw_sd',
       'POC_raw_unit', 'POC_raw_qc', 'Flux_total', 'Flux_total_unit',
       'POC_short_name_PANGAEA', 'POC_name_PANGAEA', 'ID'],
      dtype='object')
"""
#remove if it is on land
pocdf = alldf[alldf["Land_flag"] == False]

#filter to only date, position, measurement and depth and method
pocdf = pocdf[['Latitude [°]', 'Longitude [°]', 'Depth [m]',
              'Date_num', 'POC_converted [mg/m-2/d-1]','General_category']]

#drop any missing values
pocdf = pocdf.dropna(subset =['Latitude [°]', 'Longitude [°]', 'Depth [m]',
               'Date_num', 'POC_converted [mg/m-2/d-1]','General_category'] )
print(pocdf.shape)

#%%
#set the dataframe time to the right time
#"Date_num" represents the UNIX time (seconds since 1970-01-01)

#AS A QUICK FIX BC fromtimestamp doesnt except negative numbers and chla data starts at 1997:
pocdf = pocdf[pocdf["Date_num"] >= 0]#.reset_index(drop = True)
#remove this if using climatology/older data for analysis

#make dates column
dates = []
for i, r in pocdf.iterrows():
    dates.append(dt.datetime.fromtimestamp(r["Date_num"]))
pocdf["Date"] = dates

#the sst finishes in 2022 and the chla starts in 1997.
#bound the dataset to this.
pocdf = pocdf[pocdf["Date"] >= dt.datetime(1997, 9, 1)]
pocdf = pocdf[pocdf["Date"] < dt.datetime(2022, 7, 1)]#.reset_index(drop = True)
print(pocdf.shape)
#%%
#set up lists to append data to
ssts = []
chlas = []

#inport SST data for matching
sst = xr.open_dataset(f"{fp}input_data/sst/correct_lon_SST.nc")["SST"]

#slow loop for peace of mind it is working correctly
#for each row in the pocdf daatframe
for i, r in  pocdf.iterrows():
    #print(i)
    #get lat and lon position
    lat = r["Latitude [°]"]
    lon = r["Longitude [°]"]

    #set the day as the middate of the month to match to SST
    day = r["Date"]
    sday = dt.datetime(day.year, day.month, 15)

    #set the latitude and longitude to match the middle of the grid cell
    if lon == 180: #does not read in properly if not
        lon = 179
    
    #make new lat and lon
    tlat = temp_match(lat)
    tlon = temp_match(lon)

    #extract the sst value for this datapoint
    lattemp = sst.sel(lat = tlat) #select the values at the latitude
    lontemp = lattemp.sel(lon = tlon) #from this select the values at the longitude
    daytemp = lontemp.sel(time = sday, method = "nearest") #from this select the values at the date
                                        #nearest is used because the date can vary from th 14-16th
    
    #add the mean of this to the sst list
    ssts.append(np.nanmean(daytemp.values))

    #open the correct chla file
    chla = xr.open_dataset(f"{fp}input_data/oc_cci/occci_chla_{day.year}.nc")["chlor_a"]

    #set the date to match the chla days
    #chla is always the first of the month apart from the first month
    cday = dt.datetime(day.year, day.month, 1)
    chla = chla.sel(time = cday, method = "nearest")

    #find mean oc-cci for the degree cell it is in
    #select by latitude
    latchla = chla.sel(lat = slice(lat+ 0.5, lat- 0.5))
    if latchla.shape[0] == 0: #make sure it has shape, if it doesnt, do the other way
        latchla = chla.sel(lat = slice(lat- 0.5, lat+ 0.5))
    #select by longitude
    lonchla = latchla.sel(lon = slice(lon- 0.5, lon+ 0.5))
    if lonchla.shape[1] == 0:  #make sure it has shape, if it doesnt, do the other way
        lonchla = latchla.sel(lon = slice(lon+ 0.5, lon- 0.5))
    if lonchla.shape[0] * lonchla.shape[1] != 24*24:
        print(i, lonchla.shape)

    #add the mean to the dataset
    chlas.append(np.nanmean(lonchla.values))

#%%
pocdf["SST"] = ssts
pocdf["Chla"] = chlas
#select wanted columns and rename them

pocdf= pocdf.rename(columns= {"Latitude [°]":"Latitude", "Longitude [°]":"Longitude", 
                            "General_category":"Method", 'Depth [m]':"Depth", 'Date_num':"Unix_Date",
                             'POC_converted [mg/m-2/d-1]':"POC"})

#%%
#put bounds on the data
pocdf = pocdf[pocdf["Depth"] > 100]
pocdf = pocdf[pocdf["Chla"] >= 0.02]
pocdf = pocdf[pocdf["POC"] >= 0.1]
pocdf = pocdf[pocdf["SST"] > -1.79]

#subset to columns
poc = pocdf[["Latitude""Longitude","Method", "Depth", "Unix_Date", "Date",
            "POC", "SST", "Chla", 'Instrument']]

#drop any nans and no data
poc = poc.dropna()
poc = poc.drop_duplicates()

#log transform :)
poc["log_Chla"] = np.log(poc["Chla"]*1000)
poc["log_SST"] = np.log(poc["SST"] + 1.79)
poc["log_Depth"] = np.log(poc["Depth"])
poc["log_POC"] = np.log(poc["POC"])

#save
poc.to_csv(f"{fp}/merged_POC_190526.csv", index = False)

# %%
#Plot extent of the data (S1)
alldf = pd.read_csv(f"{fp}/Global_POC_Database_2026-03-10.csv")
pocdf = alldf[alldf["Land_flag"] == False]
pocdf = pocdf[['Latitude [°]', 'Longitude [°]', 'Depth [m]',
                'Date_num', 'POC_converted [mg/m-2/d-1]','General_category']]
pocdf = pocdf.dropna()
poc = pd.read_csv(f"{fp}merged_POC_190526.csv")

fig = plt.figure(figsize=(8, 6), dpi = 300)
ax = plt.axes(projection=ccrs.PlateCarree())
land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#D2D2D2") 
ax.add_feature(land, zorder=0)
sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                   edgecolor='face',facecolor="#99C2E8")
ax.add_feature(sea, zorder=0)
ax.coastlines(resolution='110m', alpha=0.3)
ax.gridlines(draw_labels = True, linestyle = "--", color = "#454545", alpha  = 0.5)
#add points
ax.scatter(pocdf["Longitude [°]"], pocdf["Latitude [°]"],transform=ccrs.PlateCarree(), 
            s=15, color="#000000", linewidth = 0.2, label = f"Available Data n = {alldf.shape[0]}")
ax.scatter(poc["Longitude"], poc["Latitude"],transform=ccrs.PlateCarree(), 
            s=12, color="#c96f00", linewidth = 0, label = f"Matched Data n = {poc.shape[0]}")
plt.legend(ncols =2, bbox_to_anchor=(0.9, -0.05))
#plt.savefig(f"{fp}figs/S1_dataextent_map.png")
plt.show()

# %%
#Split for a random 10% and 90% for model testing.
poc = pd.read_csv(f"{fp}/merged_POC_190526.csv")
for x in range(10):
    tenperc = np.random.randint(poc.shape[0], size = round((poc.shape[0]*0.1)))
    poc90 = poc.drop(index = tenperc)
    poc90.to_csv(f"{fp}percents/merged_POC_90_{x}.csv", index = False)
    poc10 = poc.iloc[tenperc]
    poc10.to_csv(f"{fp}percents/merged_POC_10_{x}.csv", index = False)
