#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from glob import glob
import LISA_Curves.LISA as li
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from natsort import natsorted
import numpy as np
import os
import statsmodels.api as sm

from amuse.ext.orbital_elements import orbital_elements_from_binary
from amuse.lab import units, constants, Particles

from mathematical_functions import *
from plotter_functions import PlotterSetup


class Plotters(object):
    """Class to plot simulation results"""
    def __init__(self):
        self.format_plot = PlotterSetup()
        self.cmap = plt.get_cmap('Blues')
        
        self.tend = 2.545100089670521 | units.kyr
        self.dt = 0.1 | units.yr
        
        # Load relevant data files
        self.log_files = natsorted(glob("outputs/*.log"))
        self.energy_files = natsorted(glob("outputs/*.energies"))
        self.env_files = natsorted(glob("outputs/*.out"))
        self.initial_files = [ ]
        for file in self.log_files:
            fname = file.split("/")[-1][:-4]+".txt"
            initial_summary = os.path.join("Merger_Parameters", fname)
            self.initial_files.append(initial_summary)  

        self.N = None
        # Merging parameters
        self.sim_steps = [[ ] for i in range(7)]  
        self.merger_ecc = [[ ] for i in range(7)]  
        self.merger_sem = [[ ] for i in range(7)]  
        self.merger_type = [[ ] for i in range(7)]
        self.merger_fno = [[ ] for i in range(7)]
        
        # All binary interactions
        self.semi_major = [[ ] for i in range(7)]  
        self.eccentricity = [[ ] for i in range(7)]
        
        # GW properties
        self.GW_freq = [[ ] for i in range(7)]
        self.GW_strain = [[ ] for i in range(7)]
        
    def extract_orbital_parameters(self, p1, p2):
        """Extract the orbital parameters of the binary system"""
        binary = Particles(2)
        
        binary[0].mass = p1[0] | units.kg
        binary[0].position = p1[1:4] | units.m
        binary[0].velocity = p1[4:7] | units.ms
        
        binary[1].mass = p2[0] | units.kg
        binary[1].position = p2[1:4] | units.m
        binary[1].velocity = p2[4:7] | units.ms
        kepler_elements = orbital_elements_from_binary(binary, G=constants.G)
        
        return kepler_elements, binary
    
    def extract_merger_parameters(self, final_dt):
        """Save orbital parameters between merging particles at final dt
        Inputs:
        final_dt:  Boolean to extract the final dt or final 100 snapshots
        """
        file_no = 0
        for data_file, output in zip(self.env_files, self.log_files):
            if os.path.getsize(output) > 0:
                file_no += 1
                with open(output, 'r') as df:
                    data_file = np.loadtxt(data_file)
                    data_values = df.readlines()
                    
                    merger_idx = data_values[-1].split("=")[-1]
                    merger_idx_1 = merger_idx[1]
                    merger_idx_2 = merger_idx[3:]
                    
                    if merger_idx_1 == "-" or merger_idx_2 == "-":
                        print(f"Curious? No mergers {output}")
                        continue
                    
                    merger_idx_1 = int(merger_idx_1)
                    merger_idx_2 = int(merger_idx_2)
                    
                    sim_pop = int(data_values[0].split("=")[1])
                    sim_dt = int(len(data_file) / sim_pop)
                    
                    if (final_dt):
                        merger1_row = (sim_dt-1)*sim_pop + merger_idx_1
                        merger2_row = (sim_dt-1)*sim_pop + merger_idx_2
                        array_idx = sim_pop // 5 - 2
                        
                        p1 = data_file[merger1_row, :]
                        p2 = data_file[merger2_row, :]
                        kepler_elements, dummy = self.extract_orbital_parameters(p1, p2)
                        
                        self.merger_sem[array_idx].append(kepler_elements[2].value_in(units.au))
                        self.merger_ecc[array_idx].append(kepler_elements[3])
                    
                    else:
                        offset = min(200, sim_dt)
                        for dt in range(offset):
                            merger1_row = (sim_dt-offset+dt)*sim_pop + merger_idx_1
                            merger2_row = (sim_dt-offset+dt)*sim_pop + merger_idx_2
                            array_idx = sim_pop // 5 - 2
                            
                            p1 = data_file[merger1_row, :]
                            p2 = data_file[merger2_row, :]
                            kepler_elements, binary = self.extract_orbital_parameters(p1, p2)
                            semimajor = kepler_elements[2]
                            eccentricity = kepler_elements[3]
                            
                            #GW calculations
                            strain_val = gw_strain(semimajor, eccentricity, binary[0].mass, binary[1].mass)
                            freq_val = gw_freq(semimajor, eccentricity, binary[0].mass, binary[1].mass).value_in(units.Hz)
                            
                            self.semi_major[array_idx].append(semimajor.value_in(units.au))
                            self.eccentricity[array_idx].append(1 - eccentricity)
                            self.GW_freq[array_idx].append(freq_val)
                            self.GW_strain[array_idx].append(strain_val)
                            self.merger_fno[array_idx].append(file_no)
                            
                        print(eccentricity)
        
    def final_ecc_CDF(self):
        self.extract_merger_parameters(final_dt=True)
        
        fig, ax = plt.subplots()
        ax.set_ylabel(r"$f_<$", fontsize=self.format_plot.font_size)
        ax.set_xlabel(r"$\log_{10}(1-e(t_{\rm coll}))$", fontsize=self.format_plot.font_size)
        self.format_plot.tickers(ax)
        
        # Plot per population
        for i, pop in enumerate(self.merger_ecc):
            ecc = np.array([1-item for item in pop])
            cdf_ecc_x, cdf_ecc_y = self.format_plot.cdf_maker(ecc)
            colour = self.cmap((i+3)/10)
            
            ax.plot(np.log10(cdf_ecc_x), cdf_ecc_y, color=colour, 
                    label=r"$N = $"+str(5*(i+2))
            )
                
        all_ecc = np.array([1-item for sublist in self.merger_ecc for item in sublist])
        cdf_ecc_x, cdf_ecc_y = self.format_plot.cdf_maker(all_ecc)
        ax.plot(np.log10(cdf_ecc_x), cdf_ecc_y, color="black", linewidth=4, label="All")
                
        ax.legend(fontsize=self.format_plot.font_size)
        ax.set_xlim(np.log10(all_ecc).min(), np.log10(all_ecc).max())
        ax.set_ylim(0,1)
        plt.savefig("plotters/final_ecc_cdf.png", dpi=300, bbox_inches='tight')
        plt.clf()
       
    def frequency_strain(self):
        """Plot the frequency vs. strain diagram of mergers"""
        self.extract_merger_parameters(final_dt=False)
        ecc_normalise = plt.Normalize(-8, 0)
        all_freq = [ ]
        all_strain = []
        
        # Set up future GW interferometers overlay
        lisa = li.LISA() 
        x_temp = np.linspace(1e-5, 1, 1000)
        Sn = lisa.Sn(x_temp)
        
        Ares = np.load(os.path.join(os.getcwd(), 'plotters/SGWBProbe/files/S_h_muAres_nofgs.npz'))
        Ares_freq = Ares['x']
        Ares_strain = Ares['y']

        fig = plt.figure(figsize=(9,7))
        gs = fig.add_gridspec(2, 2,  width_ratios=(4, 2), height_ratios=(2, 4),
                              left=0.1, right=0.9, bottom=0.1, top=0.9,
                              wspace=0.15, hspace=0.1
                              )
        ax = fig.add_subplot(gs[1, 0])
        ax1 = fig.add_subplot(gs[0, 0], sharex=ax)
        ax2 = fig.add_subplot(gs[1, 1], sharey=ax)
        ax.set_xlabel(r"$\log_{10} f$ [Hz]", fontsize=self.format_plot.font_size)
        ax.set_ylabel(r"$\log_{10} h$", fontsize=self.format_plot.font_size)
        ax1.set_ylabel(r"$\log_{10} h$", fontsize=self.format_plot.font_size)
        ax2.set_xlabel(r"$\log_{10} f$ [Hz]", fontsize=self.format_plot.font_size)
        ax.plot(np.log10(x_temp), np.log10(np.sqrt(x_temp*Sn)), color='black', linewidth=1.5)
        ax.plot(np.log10(Ares_freq), np.log10(np.sqrt(Ares_freq*Ares_strain)), color='black', linewidth=1.5)
        ax.text(-1, -20.5, r'$\mu$Ares', fontsize=self.format_plot.font_size, rotation=30, color='black')
        ax.text(-1, -19, 'LISA',fontsize=self.format_plot.font_size, rotation=30, color='black')
        
        norm = ecc_normalise
        sm_colour = plt.cm.ScalarMappable(cmap=self.cmap, norm=norm)
        sm_colour.set_array([])

        # Define color map and normalization for edge colors
        edge_colors = self.format_plot.pop_colours
        num_colors = len(edge_colors)
        edge_cmap = mcolors.ListedColormap(edge_colors)
        edge_color_boundaries = np.arange(num_colors + 1) - 0.5
        edge_color_norm = mcolors.BoundaryNorm(boundaries=edge_color_boundaries, ncolors=num_colors)
        sm_edge = plt.cm.ScalarMappable(cmap=edge_cmap, norm=edge_color_norm)
        sm_edge.set_array([])

        for ax_ in [ax, ax1, ax2]:
            self.format_plot.tickers(ax_)
            
        for i in range(7):
            freq = np.asarray(self.GW_freq[i])
            strain = np.asarray(self.GW_strain[i])
            ecc_vals = np.asarray(self.eccentricity[i])
            
            for sim in np.unique(self.merger_fno[i]):
                indices = np.where(np.array(self.merger_fno[i]) == sim)[0]
                ecc0 = np.asarray([np.log10(ecc_vals[indices[-1]]) for i in indices])
                edge_colour = edge_colors[i]
                
                scatter = ax.scatter(np.log10(freq[indices]), np.log10(strain[indices]), 
                                     c=ecc0, norm=ecc_normalise, cmap=self.cmap, s=10,
                                     )
                ax.scatter(np.log10(freq[indices][-1]), np.log10(strain[indices][-1]), 
                           c=ecc0[-1], norm=ecc_normalise, marker="X", s=100, 
                           edgecolors=edge_colour, linewidth=1.5, cmap=self.cmap
                           )
                ax.plot(np.log10(freq[indices]), np.log10(strain[indices]), color="black", alpha=0.2)
                
            all_freq = np.concatenate((all_freq, freq), axis=None)
            all_strain = np.concatenate((all_strain, strain), axis=None)
            
        kde = sm.nonparametric.KDEUnivariate(np.log10(all_freq))
        kde.fit()
        kde.density = kde.density/max(kde.density)
        ax1.plot(kde.support, kde.density, color="black")
        ax1.fill_between(kde.support, kde.density, alpha=0.4, color="black")
        
        kde = sm.nonparametric.KDEUnivariate(np.log10(all_strain))
        kde.fit()
        kde.density = kde.density/max(kde.density)
        ax2.plot(kde.density, kde.support, color="black")
        ax2.fill_between(kde.density, kde.support, alpha=0.4, color="black")
            
        # Eccentricity colours
        color_bar = plt.colorbar(sm_colour, ax=ax1, orientation='vertical')
        color_bar.set_label(r"$\log_{10}(1-e(t_{\rm coll}))$", fontsize=self.format_plot.font_size)

        # Population colours
        edge_color_bar = plt.colorbar(sm_edge, ax=ax1, orientation='vertical')
        edge_color_bar.set_label(r"$N$", fontsize=self.format_plot.font_size)        
        ticks = [0, 1, 2, 3, 4, 5, 6]  # Positions for each color category
        tick_labels = ['10', '15', '20', '25', '30', '35', '40']  # Custom labels
        edge_color_bar.set_ticks(ticks)
        edge_color_bar.set_ticklabels(tick_labels)
        
        ax.set_ylim(-26.9, -17.02)
        ax1.set_ylim(0, 1.02)
        ax2.set_ylim(-26.9, -17.02)
        ax.set_xlim(-3.98, 0.2)
        ax1.set_xlim(-3.98, 0.2)
        ax2.set_xlim(0, 1.02)
        plt.savefig("plotters/GW_strain_freq.png", dpi=300, bbox_inches='tight')
        plt.clf()
        STOP
         
    def tGW_comparison(self):
        """Plot the comparison between the analytical Peters and BrutusPN tGW"""
        Nhyperbolic = 0
        nhead = 0
        rSBH = 3*(2*constants.G*(4.001e6|units.MSun))/(constants.c**2)
        
        cluster_pop = [ ]
        sim_tGW = [ ]
        peters_tGW = [ ]
            
        for sim_data, initial_data in zip(self.log_files, self.initial_files):
            if os.path.getsize(sim_data) > 0:
                with open(initial_data, 'r') as df:
                    data_values = df.readlines()
                    prd_tGW_df = float(data_values[-1].split(":")[1][:-3])
                    prd_ecc = float(data_values[2].split(":")[1][:-1])
                    prd_sem = float(data_values[1].split(":")[1][:-3]) | units.au
                    
                    rp = prd_sem*(1-prd_ecc)
                    if not np.isinf(prd_tGW_df) and rp > rSBH:
                        prd_tGW_df *= 1 | units.yr
                        peters_tGW.append(prd_tGW_df.value_in(units.kyr))

                        with open(sim_data, 'r') as df:
                            data_values = df.readlines()
                            Ncluster = int(data_values[0].split("=")[1])
                            sim_tGW_df = float(data_values[3].split("=")[1])
                            sim_tGW_df *= self.tend
                            
                            cluster_pop.append(Ncluster)
                            sim_tGW.append(sim_tGW_df.value_in(units.kyr))
                    
                    if prd_ecc > 1:
                        Nhyperbolic += 1
                    if rp <= rSBH:
                        nhead += 1

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

        fig, ax = plt.subplots()
        self.format_plot.tickers(ax)
        ax.set_xlabel(r"$N_{\rm IMBH}$", fontsize=self.format_plot.font_size)
        ax.set_ylabel(r"$\log_{10} t_{\rm BrutusPN} / t_{\rm Peters}$", fontsize=self.format_plot.font_size)
        ax.scatter(np.unique(cluster_pop), np.log10(merger_ratio), color="black",)
        ax.scatter(np.unique(cluster_pop), np.log10(merger_lower), color="black", marker="_")
        ax.scatter(np.unique(cluster_pop), np.log10(merger_higher), color="black", marker="_")
        ax.plot([np.unique(cluster_pop), np.unique(cluster_pop)],
                [np.log10(merger_lower), np.log10(merger_higher)], color="black"
                )
        ax.axhline(0, color="black", linestyle=":")
        plt.show()

        print(f"Initial conditions had {Nhyperbolic} hyperbolic mergers")
        print(f"Initial conditions had {nhead} head-ons mergers")
        
dd = Plotters()
dd.final_ecc_CDF()
# dd.tGW_comparison()
dd.frequency_strain()