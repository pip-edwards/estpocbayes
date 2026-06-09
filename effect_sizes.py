"""
Sample code for GRL submission: Capturing the Global Variability of Marine Particulate Organic Carbon Flux: A Hierarchical Bayesian Approach
Code plots effect sizes (Figure 3)

Created by: Pippa Edwards 
"""
#%%
#Plotting the caterpillars

#Packages
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
sns.set_context("paper")

#set fp
fp = ""
os.chdir(fp)

#name file
file = []


# %%
#Plot summative effect sizes and the complicated caterpillar plot
#Transform the data by multiply by SD

#REMOVED FROM LOOP AS NOW ONLY WORKING WITH ONE FILE

#load in betas
betas = pd.read_csv(glob.glob(f"{fp}{file}/*_beta_vals.csv")[0])
#load in gammas
gammas = pd.read_csv(glob.glob(f"{fp}{file}/*_gamma_vals.csv")[0])

#rename the columns
betas["T"] = betas["beta.2"] 
betas["C"] = betas["beta.3"] 
betas["Z"] = betas["beta.4"] 
betas["T*C"] = betas["beta.5"]
betas["T*Z"] = betas["beta.6"] 
betas["C*Z"] = betas["beta.7"]
gammas["T"] = gammas["gamma.2"] 
gammas["C"] = gammas["gamma.3"] 
gammas["Z"] = gammas["gamma.4"] 
gammas["T*C"] = gammas["gamma.5"]
gammas["T*Z"] = gammas["gamma.6"] 
gammas["C*Z"] = gammas["gamma.7"]

#set up new data frames with the nicer names
bbetas = betas[["T","C","Z","T*C","T*Z","C*Z"]]
ggammas = gammas[["T","C","Z","T*C","T*Z","C*Z"]]
#make a list of column names
columns = ["T","C","Z","T*C","T*Z","C*Z"]

#load in original data
pocdf = pd.read_csv(f"{fp}merged_POC_190526.csv")
#set arrays of predictor variables 
sst = np.array(pocdf["log_SST"])
chla = np.array(pocdf["log_Chla"])
depth = np.array(pocdf["log_Depth"])

#make datasets
call =  (bbetas["C"] + np.median(np.multiply.outer(np.array(bbetas["T*C"]),sst), axis = 1) + np.median(np.multiply.outer(np.array(bbetas["C*Z"]),depth), axis = 1))*np.std(pocdf["log_Chla"])
sall = (bbetas["T"] + np.mean(np.multiply.outer(np.array(bbetas["T*C"]),chla), axis = 1) + np.mean(np.multiply.outer(np.array(bbetas["T*Z"]),depth), axis = 1))*np.std(pocdf["log_SST"])
dall = (bbetas["Z"] +  np.mean(np.multiply.outer(np.array(bbetas["C*Z"]),chla), axis = 1) +  np.mean(np.multiply.outer(np.array(bbetas["T*Z"]),sst),axis = 1))*np.std(pocdf["log_Depth"])
#summarise all betas
balls = pd.DataFrame({"T":sall,"C":call, "Z":dall})

#make datasets 
call =  (ggammas["C"] + np.median(np.multiply.outer(np.array(ggammas["T*C"]),sst), axis = 1) + np.median(np.multiply.outer(np.array(ggammas["C*Z"]),depth), axis = 1))*np.std(pocdf["log_Chla"])
sall = (ggammas["T"] + np.mean(np.multiply.outer(np.array(ggammas["T*C"]),chla), axis = 1) + np.mean(np.multiply.outer(np.array(ggammas["T*Z"]),depth), axis = 1))*np.std(pocdf["log_SST"])
dall = (ggammas["Z"] +  np.mean(np.multiply.outer(np.array(ggammas["C*Z"]),chla), axis = 1) +  np.mean(np.multiply.outer(np.array(ggammas["T*Z"]),sst),axis = 1))*np.std(pocdf["log_Depth"])
#summarise all gammas
galls = pd.DataFrame({"T":sall,"C":call, "Z":dall})

#%%
#make plots
#do in two halves for ease
cs1 = ["#ff7433","#62a87c","#3969d1"]
cs2 = ["#ffbe33","#9ad4b0","#8CC1E0","#cabf59","#beabf3","#CC89A8"]

#BETAS
yticks1 = [r"$\beta_{T}$",r"$\beta_{C}$",r"$\beta_{Z}$"]
yticks2 = [r"$\beta_{T0}$",r"$\beta_{C0}$",r"$\beta_{Z0}$",r"$\beta_{TC}$",r"$\beta_{TZ}$",r"$\beta_{CZ}$",]

fig,ax = plt.subplots(2,1, figsize = (6,7), dpi = 300)
#plot violins
sns.violinplot(data = balls, palette = cs1, inner = None, orient = "h", gap = -0.3,
            linewidth = 0.8, linecolor = "#333333", alpha = 0.8, ax = ax[0])
sns.violinplot(data = bbetas, palette = cs2, inner = None, orient = "h", gap = -0.3,
            linewidth = 0.8, linecolor = "#333333", alpha = 0.8)

#plot median lines
ax[0].scatter(np.median(balls, axis = 0),np.arange(0,3,1), c= "#333333", s = 10, marker  = "|", alpha = 0.8)
ax[1].scatter(np.median(bbetas, axis = 0),np.arange(0,6,1), c= "#333333", s = 10, marker  = "|", alpha = 0.8)

#sort axis
ax[0].set_yticks(np.arange(0,3,1), yticks1, size = 15)
ax[1].set_yticks(np.arange(0,6,1), yticks2, size = 13)

ax[0].set_xticks(np.arange(-0.8, 0.9,0.2))
ax[1].set_xlim([-1.45,1.25])
ax[1].set_xticks(np.arange(-1.5, 1.3,0.3))

for i in range(2):
    ax[i].axline((0,0), (0,1),  ls = "--", color = "black")
    ax[i].set_xlabel("Effect Size on μ", size = 12)
    ax[i].tick_params(axis='x', labelsize=10) 
#plt.show()
#plt.savefig(f"{fp}/figs/{file}/effect_sizes_beta.png")#, transparent = True)


#GAMMAS
yticks3 = [r"$\gamma_{T}$",r"$\gamma_{C}$",r"$\gamma_{Z}$"]
yticks4 = [r"$\gamma_{T0}$",r"$\gamma_{C0}$",r"$\gamma_{Z0}$",r"$\gamma_{TC}$",r"$\gamma_{TZ}$",r"$\gamma_{CZ}$",]

fig,ax = plt.subplots(2,1, figsize = (6,7), dpi = 300)
#plot violins
sns.violinplot(data = galls, palette = cs1, inner = None, orient = "h", gap = -0.3,
            linewidth = 0.8, linecolor = "#333333", alpha = 0.8, ax = ax[0])

sns.violinplot(data = ggammas, palette = cs2, inner = None, orient = "h",gap = -0.3,
            linewidth = 0.8, linecolor = "#333333", alpha = 0.8)#,  

#plot median lines
ax[0].scatter(np.median(galls, axis = 0),np.arange(0,3,1), c= "#333333", s = 10, marker  = "|", alpha = 0.8)
ax[1].scatter(np.median(ggammas, axis = 0),np.arange(0,6,1), c= "#333333", s = 10, marker  = "|", alpha = 0.8)

#sort axis
for i in range(2):
    ax[i].axline((0,0), (0,1),  ls = "--", color = "black")
    ax[i].set_xlabel("Effect Size on σ", size = 12)
    ax[i].tick_params(axis='x', labelsize=10) 

ax[0].set_yticks(np.arange(0,3,1), yticks3, size = 15)
ax[1].set_yticks(np.arange(0,6,1), yticks4, size = 13)

ax[0].set_xticks(np.arange(-0.2, 0.105,0.05))
ax[0].set_xlim([-0.2, 0.11])
ax[1].set_xlim([-1.2,0.8])
ax[1].set_xticks(np.arange(-1.2, 0.85,0.2))
#plt.show()
#plt.savefig(f"{fp}/figs/{file}/effect_sizes_gamma.png")

# %%
