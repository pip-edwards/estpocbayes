"""
Sample code for GRL submission: Capturing the Global Variability of Marine Particulate Organic Carbon Flux: A Hierarchical Bayesian Approach
Creates global estimate of POC flux
Plots figure 4 + supplementary maps.

Created by: Pippa Edwards
"""
#packages
import pandas as pd
import numpy as np
import os
import xarray as xr
from pyproj import Geod
import pandas as pd
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import os
from cmcrameri import cm
from scipy.stats import linregress
import cartopy.feature as cfeature
import glob

fp =""
file = ""
#map of where there is data (for the white line on maps)

allsum = xr.open_dataset(f"{fp}input_data/months_of_data.nc")["Months of Data"]
#%%
#import coefficient datas:
#open beta and gamma datasets 
betas = pd.read_csv(glob.glob(f"{fp}/{file}/*_beta_vals.csv")[0])
gammas = pd.read_csv(glob.glob(f"{fp}/{file}/*_gamma_vals.csv")[0])

#set up vars for each gamma (to make code lines shorter)
gamma1 = np.mean(gammas["gamma.1"])
gamma2 = np.mean(gammas["gamma.2"])
gamma3 = np.mean(gammas["gamma.3"])
gamma4 = np.mean(gammas["gamma.4"])
gamma5 = np.mean(gammas["gamma.5"])
gamma6 = np.mean(gammas["gamma.6"])
gamma7 = np.mean(gammas["gamma.7"])

beta1 = np.mean(betas["beta.1"])
beta2 = np.mean(betas["beta.2"])
beta3 = np.mean(betas["beta.3"])
beta4 = np.mean(betas["beta.4"])
beta5 = np.mean(betas["beta.5"])
beta6 = np.mean(betas["beta.6"])
beta7 = np.mean(betas["beta.7"])
#%%
####################################### FIGURES
#import depth
depth = xr.open_dataset(f"{fp}input_data/bathy/depth100_map.nc")["depth"]
depth = depth.sortby("lat")
depth = np.log(depth)

#import SST and Chl-a climatology and transform
chla_ntd = xr.open_dataset(f"{fp}input_data/occci/occci_overall_climatology.nc")["chlor_a"]
chla_ntd = np.log(chla_ntd * 1000)
sst_ntd = xr.open_dataset(f"{fp}input_data/sst/SST_overall_climatology.nc")["SST"]
sst_ntd = np.log(sst_ntd + 1.79)

#set up lats and lons for mapping purposes
lats = sst_ntd["lat"]
lons = sst_ntd["lon"]

#make distribution for mu, sigma and mean poc
mumap = beta1*np.ones((180,360)) + sst_ntd*beta2 + chla_ntd*beta3 + depth*beta4 + sst_ntd*chla_ntd*beta5 + sst_ntd*depth*beta6 + chla_ntd*depth*beta7
sigmap =  gamma1*np.ones((180,360)) + sst_ntd*gamma2 + chla_ntd*gamma3 + depth*gamma4 + sst_ntd*chla_ntd*gamma5 + sst_ntd*depth*gamma6 + chla_ntd*depth*gamma7
pocmap = mumap + (sigmap**2)/2
#%
fig = plt.figure(figsize=(10, 6), dpi = 300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines('110m', alpha=0.1)
ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")

#MEAN POC map
pocrange = np.arange(2.8, 6, 0.4)
contour = ax.contourf(lons, lats, pocmap, transform=ccrs.PlateCarree(),
                    cmap= cm.lipari,extend = "both", levels = pocrange)

cbar = plt.colorbar(contour, ax=ax, orientation='vertical', 
                    shrink=0.8, pad=0.05, extend = "both")
cbar.set_label(f"Mean ln(POC) (mg/m²/day)", size = 14 )#, weight = "bold")
land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#555555") 
contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                    linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
ax.add_feature(land, zorder=2)
sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                edgecolor='face',facecolor="#ADADAD")
ax.add_feature(sea, zorder=0)
cbar.ax.tick_params(labelsize=14)
plt.savefig(f"{fp}/figs/{file}/climatology_poc_map.png", transparent = True)
plt.show()

#%

#MU MAP
murange = np.arange(1.8, 5, 0.4)
fig = plt.figure(figsize=(10, 6), dpi = 300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines('110m', alpha=0.1)
ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
contour = ax.contourf(lons, lats, mumap, transform=ccrs.PlateCarree(),
                    cmap= cm.lipari,extend = "both", levels = murange)
cbar = plt.colorbar(contour, ax=ax, orientation='vertical', 
                    shrink=0.8, pad=0.05, extend = "both")
cbar.set_label(f"μ", size = 14)#, weight = "bold")
land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#555555") 
contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                    linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
ax.add_feature(land, zorder=2)
sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                edgecolor='face',facecolor="#ADADAD")
ax.add_feature(sea, zorder=0)
cbar.ax.tick_params(labelsize=14)
plt.savefig(f"{fp}/figs/{file}/climatology_mu_map.png", transparent = True)
plt.show()
#%
#SIG map
fig = plt.figure(figsize=(10, 6), dpi = 300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines('110m', alpha=0.1)
ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
sigrange = np.arange(1.1,2.1,0.1)
contour = ax.contourf(lons, lats, sigmap, transform=ccrs.PlateCarree(),
                    cmap= cm.lipari,extend = "both", levels = sigrange)
cbar = plt.colorbar(contour, ax=ax, orientation='vertical', 
                    shrink=0.8, pad=0.05, extend = "both")
cbar.set_label(f"σ", size = 14)#, weight = "bold")
land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#555555") 
contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                    linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
ax.add_feature(land, zorder=2)
sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                edgecolor='face',facecolor="#ADADAD")
ax.add_feature(sea, zorder=0)
cbar.ax.tick_params(labelsize=14)
plt.savefig(f"{fp}/figs/{file}/climatology_sig_map.png", transparent = True)
plt.show()

#%
#BetaZ and GammaZ term generation
betaz = beta4 + sst_ntd*beta6 + chla_ntd*beta7
gammaz =  gamma4 + sst_ntd*gamma6 + chla_ntd*gamma7

#Bz Map
fig = plt.figure(figsize=(10, 6), dpi = 300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines('110m', alpha=0.1)
ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
crange = np.arange(-1, 0.1, 0.1)
contour = ax.contourf(lons, lats, betaz, transform=ccrs.PlateCarree(),
                    cmap= cm.lipari_r,extend = "min", levels = crange)
cbar = plt.colorbar(contour, ax=ax, orientation='vertical', 
                    shrink=0.8, pad=0.05, extend = "min")
land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#555555") 
contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                    linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
ax.add_feature(land, zorder=2)
sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                edgecolor='face',facecolor="#ADADAD")
ax.add_feature(sea, zorder=0)
cbar.set_label(r"$β_z$", size = 14)
cbar.ax.tick_params(labelsize=14)
plt.savefig(f"{fp}/figs/{file}/betaz_map.png", transparent = True)
plt.show()

#Yz Map
fig = plt.figure(figsize=(10, 6), dpi = 300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines('110m', alpha=0.1)
ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
crange = np.arange(-0.3, 0.05, 0.05)
contour = ax.contourf(lons, lats, gammaz, transform=ccrs.PlateCarree(),
                    cmap= cm.lipari_r,extend = "both", levels = crange)
cbar = plt.colorbar(contour, ax=ax, orientation='vertical', 
                    shrink=0.8, pad=0.05, extend = "min")
land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#555555") 
contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                    linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
ax.add_feature(land, zorder=2)
sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                edgecolor='face',facecolor="#ADADAD")
ax.add_feature(sea, zorder=0)
cbar.set_label(r'$γ_z$', size = 14)
cbar.ax.tick_params(labelsize=14)
plt.savefig(f"{fp}/figs/{file}/gammaz_map.png", transparent = True)
plt.show()

#b map
for x in np.arange(100,1600,100):
    #print(x)

    #open depth map 
    depth = xr.open_dataset(f"{fp}input_data/bathy/depth{x}_map.nc")["depth"]
    depth = depth.sortby("lat")
    depth = np.log(depth)

    #make new calc for mean poc
    mumap = beta1*np.ones((180,360)) + sst_ntd*beta2 + chla_ntd*beta3 + depth*beta4 + sst_ntd*chla_ntd*beta5 + sst_ntd*depth*beta6 + chla_ntd*depth*beta7
    sigmap =  gamma1*np.ones((180,360)) + sst_ntd*gamma2 + chla_ntd*gamma3 + depth*gamma4 + sst_ntd*chla_ntd*gamma5 + sst_ntd*depth*gamma6 + chla_ntd*depth*gamma7
    pocmap = mumap + (sigmap**2)/2

    #add to a dataframe
    if x == 100:
        depthpoc = pocmap.expand_dims(depth = [x])
    else:
        depthpoc1 = pocmap.expand_dims(depth = [x]) #make a new one for every point that is not the first
        depthpoc = xr.concat([depthpoc, depthpoc1], dim = "depth") #add to the first one

#set up empty data array
bmap = np.ones((180,360))
#set up list of log depths for the linear regression
depths = np.log(np.arange(100,1600,100))

#for each lat and lon
for lat in range(180):
    for lon in range(360):
        #get the pocs value at each depth
        pocvals = depthpoc[:,lat,lon]
        pocvals = np.array(pocvals.values)
        #make sure it has data within it
        i = len(pocvals) - np.isnan(pocvals).sum()
        #if it doesnt, put a nan in
        if i == 0:
            bmap[lat,lon] = np.nan
        #otherwise linear regress for b
        else:
            b,_,_,_,_ = linregress(depths[:i],pocvals[:i])
            bmap[lat,lon] = b

#plot b map
fig = plt.figure(figsize=(10, 6), dpi = 300)
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.coastlines('110m', alpha=0.1)
ax.gridlines(draw_labels = True, linestyle = "--", color = "#B3B3B3")
#bmap = bmap.where(bmap<0, np.nan)
crange = np.arange(-1, -0.3, 0.1)
contour = ax.contourf(lons, lats, bmap, transform=ccrs.PlateCarree(),
                    cmap= cm.lipari_r,extend = "both", levels = crange)
cbar = plt.colorbar(contour, ax=ax, orientation='vertical', 
                    shrink=0.8, pad=0.05, extend = "both")
land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#555555") 
sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                edgecolor='face',facecolor="#ADADAD")
ax.add_feature(sea, zorder=0)
contour2 = ax.contour(lons, lats, allsum, transform=ccrs.PlateCarree(),
                    linewidths=1, cmap= cm.grayC_r, levels = np.arange(11,13,1))
ax.add_feature(land, zorder=2)
cbar.ax.tick_params(labelsize=14)
cbar.set_label(r'b', size = 14)
plt.savefig(f"{fp}review_run1905/figs/{file}/bmap_pos.png", transparent = True)
plt.show()


#%
###MONTHLY MAPS (SUPPLEMENTARY FIGURES)

#month list
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", 
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

#re-import depth
depth = xr.open_dataset(f"{fp}input_data/bathy/depth100_map.nc")["depth"]
depth = depth.sortby("lat")
depth = np.log(depth)

#Mean POC
#set up range
pocrange = np.arange(2.8, 6, 0.4)
fig, axes =  plt.subplots(nrows=3, ncols=4, figsize=(12, 6),
                        subplot_kw={'projection': ccrs.Robinson()}, dpi = 300)
fig.subplots_adjust(hspace=0.05, wspace=0.05)
for i, ax in enumerate(axes.flat):
    month = months[i]

    m = i+1
    #import chla
    chla_ntd = xr.open_dataset(f"{fp}input_data/occci/{m}_occci_monthly_regrid.nc")["chlor_a"]
    chla_ntd = np.log(chla_ntd * 1000)

    #import sst
    sst_ntd = xr.open_dataset(f"{fp}input_data/sst/SST_{m}_monthly_climatology.nc")["SST"]
    sst_ntd = np.log((sst_ntd + float(1.79)))

    mumap = beta1*np.ones((180,360)) + sst_ntd*beta2 + chla_ntd*beta3 + depth*beta4 + sst_ntd*chla_ntd*beta5 + sst_ntd*depth*beta6 + chla_ntd*depth*beta7
    sigmap =  gamma1*np.ones((180,360)) + sst_ntd*gamma2 + chla_ntd*gamma3 + depth*gamma4 + sst_ntd*chla_ntd*gamma5 + sst_ntd*depth*gamma6 + chla_ntd*depth*gamma7
    pocmap = mumap + (sigmap**2)/2

    contour = ax.contourf(lons, lats, pocmap, transform=ccrs.PlateCarree(),
                            cmap=cm.lipari, extend = "both", levels = pocrange)
    ax.set_title(month, weight = "bold", size = 13)

    ax.coastlines('110m', alpha=0.1)
    land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#555555") 
    ax.add_feature(land, zorder=1)
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                    edgecolor='face',facecolor="#ADADAD")
    ax.add_feature(sea, zorder=0)

cbar_ax = fig.add_axes([0.2, 0.08, 0.6, 0.03])#0.1, 0.08, 0.7, 0.05])  # [left, bottom, width, height]

cbar = fig.colorbar(contour, cax=cbar_ax, orientation='horizontal', extend = "both")
cbar.set_label('ln(POC) (mg/m²/day)', size = 13)
cbar.ax.tick_params(labelsize=12)
plt.savefig(f"{fp}/figs/{file}/poc_map_allmonths.png", transparent = True)


#MU
murange = np.arange(1.8, 5, 0.4)
fig, axes =  plt.subplots(nrows=3, ncols=4, figsize=(12, 6),
                        subplot_kw={'projection': ccrs.Robinson()}, dpi = 300)
fig.subplots_adjust(hspace=0.05, wspace=0.05)
for i, ax in enumerate(axes.flat):
    month = months[i]
    m = i+1
    #imoprt chla
    chla_ntd = xr.open_dataset(f"{fp}input_data/occci/{m}_occci_monthly_regrid.nc")["chlor_a"]
    chla_ntd = np.log(chla_ntd * 1000)

    #import sst
    sst_ntd = xr.open_dataset(f"{fp}input_data/sst/SST_{m}_monthly_climatology.nc")["SST"]
    sst_ntd = np.log((sst_ntd + float(1.79)))
    mumap = beta1*np.ones((180,360)) + sst_ntd*beta2 + chla_ntd*beta3 + depth*beta4 + sst_ntd*chla_ntd*beta5 + sst_ntd*depth*beta6 + chla_ntd*depth*beta7

    contour = ax.contourf(lons, lats, mumap, transform=ccrs.PlateCarree(),
                            cmap=cm.lipari, extend = "both", levels = murange)
    ax.set_title(month, weight = "bold", size = 13)

    ax.coastlines('110m', alpha=0.1)
    land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#555555") 
    ax.add_feature(land, zorder=1)
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                    edgecolor='face',facecolor="#ADADAD")
    ax.add_feature(sea, zorder=0)

cbar_ax = fig.add_axes([0.2, 0.08, 0.6, 0.03])#0.2, 0.08, 0.6, 0.03])  # [left, bottom, width, height]
cbar = fig.colorbar(contour, cax=cbar_ax, orientation='horizontal', extend = "both")
cbar.set_label('μ', size = 13)
cbar.ax.tick_params(labelsize=12)
plt.savefig(f"{fp}/figs/{file}/mu_map_allmonths.png", transparent = True)

#%
#SIGMA
sigrange = np.arange(1.1,2.1,0.1)
fig, axes =  plt.subplots(nrows=3, ncols=4, figsize=(12, 6),
                        subplot_kw={'projection': ccrs.Robinson()}, dpi = 300)
fig.subplots_adjust(hspace=0.05, wspace=0.05)
for i, ax in enumerate(axes.flat):
    month = months[i]
    m = i+1

    #imoprt chla
    chla_ntd = xr.open_dataset(f"{fp}input_data/occci/{m}_occci_monthly_regrid.nc")["chlor_a"]
    #transform chla
    chla_ntd = np.log(chla_ntd * 1000)


    #import sst
    sst_ntd = xr.open_dataset(f"{fp}input_data/sst/SST_{m}_monthly_climatology.nc")["SST"]
    sst_ntd = np.log((sst_ntd + float(1.79)))

    sigmap =  gamma1*np.ones((180,360)) + sst_ntd*gamma2 + chla_ntd*gamma3 + depth*gamma4 + sst_ntd*chla_ntd*gamma5 + sst_ntd*depth*gamma6 + chla_ntd*depth*gamma7

    contour = ax.contourf(lons, lats, sigmap, transform=ccrs.PlateCarree(),
                            cmap=cm.lipari, extend = "both", levels = sigrange)
    ax.set_title(month, weight = "bold", size = 13)

    ax.coastlines('110m', alpha=0.1)
    land = cfeature.NaturalEarthFeature('physical', 'land', '110m'
                                    ,edgecolor="#1B1B1B",facecolor="#555555") 
    ax.add_feature(land, zorder=1)
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '110m',
                                    edgecolor='face',facecolor="#ADADAD")
    ax.add_feature(sea, zorder=0)

cbar_ax = fig.add_axes([0.2, 0.08, 0.6, 0.03])#0.1, 0.08, 0.7, 0.05])  # [left, bottom, width, height]
cbar = fig.colorbar(contour, cax=cbar_ax, orientation='horizontal', extend = "both")
cbar.set_label('σ', size = 13)
cbar.ax.tick_params(labelsize=12)
plt.savefig(f"{fp}/figs/{file}/sig_map_allmonths.png", transparent = True)
#%%
###################################### ESTIMATES
#Calculate the area of each 1o cell

#import sst for lat and lon
sst_ntd = xr.open_dataset(f"{fp}input_data/sst/SST_overall_climatology.nc")["SST"]
#get lat and lon for plotting/mapping purposes
lats = sst_ntd["lat"]
lons = sst_ntd["lon"]

#select geoid
g = Geod(ellps = "WGS84") #this is the standard/most accurate eath model

#set up an empty cell area array
cell_area = np.zeros(((len(lats)), len(lons)))
for i in range(len(lats)):
    lat = lats[i]
    #print(lat)
    for j in range(len(lons)):
        lon = lons[j]
        latcell = [lat-0.5, lat-0.5, lat+0.5, lat+0.5, lat-0.5]
        loncell = [lon-0.5, lon+0.5, lon+0.5, lon-0.5, lon-0.5]  
        #fill with area
        area, _  = g.polygon_area_perimeter(loncell,latcell)
        cell_area[i,j] = area

#%%
#  QUICK GLOBAL ESTIMATE

#This is a monthly weighted estimate built off the mean beta and gamma
#set up how much each month is a fraction of the year
def get_frac(month):
    if month == 2:
        return 28/365
    elif month in [4,6,9,11]:
       return 30/365
    else:
        return 31/365

#set up a sum value to add to
allpocs = 0

#for each month
for m in range(1,13):
    #import chla
    chla_ntd = xr.open_dataset(f"{fp}input_data/occci/{m}_occci_monthly_regrid.nc")["chlor_a"]
    #transform chla to ug
    chla_ntd = np.log(chla_ntd * 1000)

    #import sst
    sst_ntd = xr.open_dataset(f"{fp}input_data/sst/SST_{m}_monthly_climatology.nc")["SST"]
    sst_ntd = np.log((sst_ntd + float(1.79)))

    #use linear equation to get mu and sigma
    mumap = beta1*np.ones((180,360)) + sst_ntd*beta2 + chla_ntd*beta3 + depth*beta4 + sst_ntd*chla_ntd*beta5 + sst_ntd*depth*beta6 + chla_ntd*depth*beta7
    sigmap =  gamma1*np.ones((180,360)) + sst_ntd*gamma2 + chla_ntd*gamma3 + depth*gamma4 + sst_ntd*chla_ntd*gamma5 + sst_ntd*depth*gamma6 + chla_ntd*depth*gamma7
    #calculate mean
    pocmap = mumap + (sigmap**2)/2

    #get fraction for the month of the total
    month_frac = get_frac(m)
    pocs = np.exp(pocmap)
    #times by cell area in m2 (mgC/m^2)
    areapocs = pocs * cell_area
    #times for 365 it is there to make yearly (mgC/m2/day)
    aypocs = areapocs*365
    #sum for the globe
    pocsum = np.sum(aypocs)
    #to gC
    pocsum = pocsum/1000
    #to PgC
    pocsum = pocsum/(10**15)
    #multiply by the fraction of the year there is data
    allpocs = allpocs + (pocsum*month_frac)
#print value
print(np.round(allpocs.values,2))
#%%
#REAL GLOBAL ESTIMATE

#This is the ugly loop to get the histograms of overall POC flux from each beta and agamma calculated.
#:D
#I could definitley optimise this
#but optimising is a step i haven't gotten to yet.

#REMEMBER ONLY RUN THIS OVERNIGHT IT WILL KILL YOUR LAPTOP IT IS SO SLOW WHY DIDNT YOU OPTIMISE THIS!!!!!!!!

#make global estimates with monthly weighting


#run through each file for each month
#define file 
f = ""
#open beta and gamma datasets 
betas = pd.read_csv(glob.glob(f"{fp}/{f}/*_beta_vals.csv")[0])
gammas = pd.read_csv(glob.glob(f"{fp}/{f}/*_gamma_vals.csv")[0])

#set up a sum value to add to (these need to be arrays because im getting a mean for each one)
pocalls = np.zeros(8000)
mualls =np.zeros(8000)
sigexpalls = np.zeros(8000)
sigalls = np.zeros(8000)

#for each month
for m in range(1,13):
    #get the month fraction
    month_frac = get_frac(m)

    #import chla
    chla_ntd = xr.open_dataset(f"{fp}input_data/occci/{m}_occci_monthly_regrid.nc")["chlor_a"]
    #transform chla depending on dataset
    chla_ntd = np.log(chla_ntd * 1000)

    #import sst
    sst_ntd = xr.open_dataset(f"{fp}input_data/sst/SST_{m}_monthly_climatology.nc")["SST"]
    sst_ntd = np.log((sst_ntd + float(1.79)))

    #(depth doesnt change)

    #set up empty dataframes for each month to add 8000 datapoints to
    allpocs = []
    allmus = []
    allsigs = []
    #cycle through each 8000 rows for beta and gamma
    for n in range(8000):
    #print(n)
    
        #select the row 
        b = betas.loc[n]
        g = gammas.loc[n]
        #use this to build estimate of mu, sig and mean POC flux at 100m
        mumap = b[0]*np.ones((180,360)) + sst_ntd*b[1] + chla_ntd*b[2] + depth*b[3] + sst_ntd*chla_ntd*b[4] +sst_ntd*depth*b[5]+ chla_ntd*depth*b[6] 
        sigmap =  g[0]*np.ones((180,360)) + sst_ntd*g[1] + chla_ntd*g[2] + depth*g[3] + sst_ntd*chla_ntd*g[4] +sst_ntd*depth*g[5]+ chla_ntd*depth*g[6] 
        pocmap = mumap + (sigmap**2)/2

        #make a globally integrated mu estimate
        #take out of logspace
        mus = np.exp(mumap)
        #times by cell area in m2
        areamus = mus * cell_area
        #make yearly
        aymus = areamus*365
        #sum for globe
        musum = np.sum(aymus)
        #to grams
        musum = musum/1000
        #to petegrams
        musum = musum/(10**15)
        #add to list
        allmus.append(musum)

        #make a globally integrated sigma estimate
        #take out of logspace
        sigs = np.exp(sigmap)
        #times by cell area in m2
        areasigs = sigs * cell_area
        #make yearly
        aysigs = areasigs*365
        #sum for globe
        sigsum = np.sum(aysigs)
        #to grams
        sigsum = sigsum/1000
        #to petagrams
        sigsum = sigsum/(10**15)
        #append to list
        allsigs.append(sigsum)

        #make a globally integrated mean estimate
        #take out of logspace
        pocs = np.exp(pocmap)
        #times by cell area in m2 (mgC/m^2)
        areapocs = pocs * cell_area
        #times for 365 it is there to make yearly (mgC/m2/day)
        aypocs = areapocs*365
        #sum for the globe
        pocsum = np.sum(aypocs)
        #to gC
        pocsum = pocsum/1000
        #to PgC
        pocsum = pocsum/(10**15)
        #append to list of summed pocs
        allpocs.append(pocsum)

    #make each list an array and times by monthly fraction
    allpocs = np.array(allpocs) *month_frac
    allmus = np.array(allmus) *month_frac
    allsigs = np.array(allsigs) *month_frac

    #add the monthly fraction to the global fraction
    pocalls = pocalls + (allpocs)
    mualls = mualls + (allmus)
    sigalls = sigalls + (allsigs)

#make a plot to visually represent this and save it 
fig, ax = plt.subplots(3,1)
fig.suptitle(f)
fig.tight_layout(h_pad = 1.2)
ax[0].hist(pocalls)
ax[0].set_title(f"poc, mean = {np.round(np.mean(pocalls),3)}, sd={np.round(np.std(pocalls),3)}")

ax[1].hist(mualls)
ax[1].set_title(f"mu, mean = {np.round(np.mean(mualls),3)}, sd={np.round(np.std(mualls),3)}")

ax[2].hist(sigalls)
ax[2].set_title(f"sig, mean = {np.round(np.mean(sigalls),3)}, sd={np.round(np.std(sigalls),3)}")
if not os.path.isdir(f"{fp}/figs/{f}"):
        os.makedirs(f"{fp}/figs/{f}")
plt.savefig(f"{fp}/figs/{f}/global_poc_est.png")
