"""
Sample code for GRL submission: Capturing the Global Variability of Marine Particulate Organic Carbon Flux: A Hierarchical Bayesian Approach
Makes the comparison table fo
Makes the plots for Figure 2 of the paper.
Plot Z-score of recreated data

Created by: Pippa Edwards
"""
#%%

#Packages
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
from scipy.stats import ks_2samp
from scipy.stats import pearsonr
from sklearn.metrics import r2_score
from cmcrameri import cm
import matplotlib.colors as colors
lipari_90 = colors.LinearSegmentedColormap.from_list(
    "lipari_90", 
    cm.lipari_r(np.linspace(0.1, 1.0, 256)))
#colour scheme :)
cs = ["#e99871", "#416c99"]
sns.set_context("paper")

#filepaths
fp = ""
os.chdir(fp)

files = ""

#%%
#make an empty dataframe to save the stats into 
comparison_stats = pd.DataFrame({"File":[],
                                 "r^2 sigma":[],
                                 "r2 no sigma":[],
                                 "r sigma":[],
                                 "r pval sigma":[],
                                 "r no sigma":[],
                                 "r pval no sigma":[],
                                 "ks sigma":[],
                                 "ks no sigma":[]})

for f in files:
    print(f)
    #split to get input data
    parts = f.split("_")
    met = parts[0]
    chl = parts[1]

    pocdf = pd.read_csv(f"{fp}merged_POC_{met}_{chl}_190526.csv")

    #### HAsH OUT THE FOLLOWING IF ALREADY GENERATED
    #(this is to save time in rerunning so mu and sigma don't have to be loaded each time)
    #get input data 
    # sigmas = pd.read_csv(glob.glob(f"{fp}{f}/*_sigma_vals.csv")[0])
    # mus = pd.read_csv(glob.glob(f"{fp}{f}/*_mu_vals.csv")[0])
    # #SET UP GENERATED LOG POC DATASETS
    # mmus = mus.mean(axis = 0).to_list()
    # msigmas = sigmas.mean(axis = 0).to_list()
    # ygen = [] #with sigma
    # ygen0 = [] #without sigma
    # for x in range(len(mmus)):
    #     y = np.random.normal(mmus[x], msigmas[x])
    #     ygen.append(y)
    #     y0 = np.random.normal(mmus[x], 0)
    #     ygen0.append(y0)
    # pocdf["ygen"] = ygen
    # pocdf["ygen0"] = ygen0
    # print(np.min(ygen), np.max(ygen))
    # pocdf.to_csv(f"{fp}{f}/{f}_input_withygen.csv")

    #load in the data 
    pocdf = pd.read_csv(f"{fp}{f}/{f}_input_withygen.csv")
    ygen = pocdf["ygen"]
    ygen0 = pocdf["ygen0"]

    ks,p1 = ks_2samp(pocdf["log_POC"],ygen)
    ks0,p0 = ks_2samp(pocdf["log_POC"],ygen0)
    #ax[0].title(f"KS with sigma = {round(ks,3)}, KS without sigma = {round(ks0,3)}")

    #make a folder if not already exisiting and save file
    if not os.path.isdir(f"{fp}/figs/{f}"):
        os.makedirs(f"{fp}/figs/{f}")
    plt.savefig(f"{fp}/figs/{f}/POC_histograms")
    plt.show()

    #print r2, ks and pcc
    print(f)
    print("r2 yegn ygen0", r2_score(pocdf["log_POC"], pocdf["ygen"]), r2_score(pocdf["log_POC"], pocdf["ygen0"])) 
    print("pcc: ygen", pearsonr(pocdf["log_POC"], pocdf["ygen"]), "ygen0",pearsonr(pocdf["log_POC"], pocdf["ygen0"]))
    print("ks ygen ygen0", ks, ks0)

    r0, p0 = pearsonr(pocdf["log_POC"], pocdf["ygen0"])
    r1, p1 = pearsonr(pocdf["log_POC"], pocdf["ygen"])

    comparison_stats.loc[len(comparison_stats)] = [f,
                                                   r2_score(pocdf["log_POC"], pocdf["ygen"]), 
                                                   r2_score(pocdf["log_POC"], pocdf["ygen0"]),
                                                   r1, p1, r0, p0, ks, ks0]


    #PLOT the figure of hisotgrams + desinty scatter
    fig, ax = plt.subplots(2,2,figsize = (8,8), sharey= "row", sharex = "row", dpi = 300,)
    plt.tight_layout(w_pad = 0.1, h_pad= 1.5)

    #histograms
    #without sigma
    sns.histplot(x = pocdf["log_POC"], alpha = 0.5, fill = True, linewidth  = 1.5,kde = True,
                color = "black", edgecolor ="black",label = "Observed POC Flux",  binwidth =0.2, ax = ax[0,1])
    sns.histplot(x = ygen0, alpha = 0.5, fill = True, edgecolor = cs[1], linewidth  = 1.5, kde = True,
                color = cs[1], label = "Modelled POC, no σ", binwidth =0.2, ax = ax[0,1])
    #with sigma
    sns.histplot(x = pocdf["log_POC"], alpha = 0.5, fill = True, linewidth  = 1.5,kde = True,
                color = "black", edgecolor ="black",label = "Observed POC Flux",  binwidth =0.2, ax = ax[0,0])
    sns.histplot(x = ygen, alpha = 0.5, fill = True, edgecolor = cs[0], linewidth  = 1.5,kde = True,
                color = cs[0], label = "Modelled POC with σ", binwidth =0.2, ax = ax[0,0])

    #hex plots
    hb = ax[1,0].hexbin(x = pocdf["log_POC"], y = pocdf["ygen"], cmap = lipari_90,
                    mincnt = 1, vmin = 0, vmax = 26, linewidths = 0.1)
    hb0 = ax[1,1].hexbin(x = pocdf["log_POC"], y = pocdf["ygen0"], cmap = lipari_90,
                    mincnt = 1, vmin = 0, vmax = 26 , linewidths = 0.1)

    #colourbar
    #cbar = fig.colorbar(hb, ax=ax, label='Density', cmap = cm.lipari)
    #cbar.ax.tick_params(labelsize=30)

    #set axis ticks/labels
    for i in range(2):
        ax[0,i].set_xlabel("ln(POC Flux) (mg/m²/day)",size = 12)
        ax[0,i].set_xticks(np.arange(-4, 9,1))
        ax[0,i].set_xlim([-4.5, 8])
        ax[0,i].legend(loc = "upper right")
        ax[0,i].tick_params(labelsize = 11)
    ax[0,0].set_ylabel("Count", size = 12)
    for i in range(2):
        ax[1,i].set_xlabel("Modelled ln(POC Flux) (mg/m²/day", size = 12)
        ax[1,i].set_xticks(np.arange(-2,9,1))
        ax[1,i].set_yticks(np.arange(-4, 9,1))
        ax[1,i].set_ylim([-4.5,8])
        ax[1,i].tick_params(labelsize = 11)
        ax[1,i].axline((0,0),(1,1), c = "black", ls = "--", linewidth = 1)
    ax[1,0].set_ylabel("Modelled ln(POC Flux) (mg/m²/day)", size = 12)

    #plt.savefig(f"{fp}/figs/{f}/POC_histoscatter.png")

comparison_stats.to_csv("histogram_comparison_stats_runs.csv", index = False)
#%%

#Z-Score
