#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from glob import glob
import matplotlib.pyplot as plt
from natsort import natsorted
import numpy as np
import os

from amuse.lab import units

dt = 2.545100089670521 | units.kyr

simulated_files = natsorted(glob("outputs/*.log"))

cluster_pop = [ ]
sim_tGW = [ ]
peters_tGW = [ ]

for exp_data in simulated_files:
    if os.path.getsize(exp_data) > 0:
        fname = exp_data.split("/")[-1][:-4]+".txt"
        out_data = os.path.join("Merger_Parameters", fname)
        with open(out_data, 'r') as df:
            data_values = df.readlines()
            prd_tGW_df = float(data_values[-1].split(":")[1][:-3])
            if not np.isinf(prd_tGW_df):
                prd_tGW_df *= 1 | units.yr
                peters_tGW.append(prd_tGW_df.value_in(units.kyr))

                with open(exp_data, 'r') as df:
                    data_values = df.readlines()
                    Ncluster = int(data_values[0].split("=")[1])
                    sim_tGW_df = float(data_values[3].split("=")[1])
                    print(sim_tGW_df*dt)
                    sim_tGW_df *= dt
                    
                    cluster_pop.append(Ncluster)
                    sim_tGW.append(sim_tGW_df.value_in(units.kyr))
        

x = np.linspace(0, max(peters_tGW), 50)
plt.plot(x, x, linestyle=":", color="black")
plt.scatter(peters_tGW, sim_tGW, c=cluster_pop)
plt.xlabel(r"$t_{\rm GW, Peters}$", fontsize=20)
plt.ylabel(r"$t_{\rm GW, BrutusPN}$", fontsize=20)
plt.ylim(0, 1.05*max(sim_tGW))
plt.xlim(0, 5)
plt.show()

merger_ratio = [ ]
merger_lower = [ ]
merger_higher = [ ]
for k in np.unique(cluster_pop):
    indices = np.where(np.array(cluster_pop) == k)
    sim_tGW_df = np.array(sim_tGW)[indices]
    pet_tGW_df = np.array(peters_tGW)[indices]
    
    ratio = sim_tGW_df/pet_tGW_df
    medL = np.percentile(ratio, 25)
    medH = np.percentile(ratio, 75)
    
    merger_ratio.append(np.median(ratio))
    merger_lower.append(medL)
    merger_higher.append(medH)

print(merger_ratio)
print(merger_lower)
print(merger_higher)

plt.scatter(np.unique(cluster_pop), np.log10(merger_ratio), color="black",)
plt.scatter(np.unique(cluster_pop), np.log10(merger_lower), color="black", marker="_")
plt.scatter(np.unique(cluster_pop), np.log10(merger_higher), color="black", marker="_")
plt.plot([np.unique(cluster_pop), np.unique(cluster_pop)],
         [np.log10(merger_lower), np.log10(merger_higher)], color="black"
         )
plt.axhline(0, color="black", linestyle=":")
plt.xlabel(r"$N_{\rm IMBH}$", fontsize=20)
plt.ylabel(r"$\log_{10} t_{\rm BrutusPN} / t_{\rm Peters}$", fontsize=20)
plt.show()