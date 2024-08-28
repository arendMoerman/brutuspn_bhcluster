#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as pt
import numpy as np
import argparse

def PlotOutput():
    parser = argparse.ArgumentParser(description="Options for BrutusPN quicklook")
    parser.add_argument("-n", "--name", type=str, help="Name of simulation result to plot.")
    parser.add_argument("-p", "--projection", type=int, default=0, help="Projection for plot. Choose 0 for xy, 1 for xz and 2 for yz. Default is 0.")
    args = parser.parse_args()

    Arr = np.loadtxt(f"./outputs/{args.name}.out")
    enr = np.loadtxt(f"./outputs/{args.name}.energies")
    
    N = int(Arr[0,-1])
    T = int(len(Arr) / N)
    m = []
    mtot = 0
    R = np.empty((N, 3, T))
    
    for i in range(N):
        m.append(Arr[i,0])

    Arr = np.delete(Arr, 0, 1)
    Arr = np.delete(Arr, -1, 1)
    
    idx_x = 0
    idx_y = 1

    if args.projection == 1:
        idx_y = 2
    elif args.projection == 2:
        idx_x = 1
        idx_y = 2

    for i in range(N):
        for k in range(T):

            R[i, 0, k] = Arr[k*N+i,0]
            R[i, 1, k] = Arr[k*N+i,1]
            R[i, 2, k] = Arr[k*N+i,2]
    
    fig, ax = pt.subplots(1, 2, width_ratios=[1, 1], height_ratios=[1])
    for i in range(N):    
        ax[0].plot(R[i,idx_x],R[i,idx_y], linewidth=1, label=('m={}'.format(round(m[i]))))
    
    for i in range(N):    
        ax[0].scatter(R[i,idx_x,-1],R[i,idx_y,-1])

    coord_labels = ["x", "y", "z"]

    ax[0].set_xlabel(f'{coord_labels[idx_x]} [m]')
    ax[0].set_ylabel(f'{coord_labels[idx_y]} [m]')
    ax[0].set_aspect("equal")
    ax[0].tick_params(which="both", direction="in")
    #ax.legend(frameon=False, prop={'size': 10})
    #pt.minorticks_on()
    #pt.savefig(fname="PlotterOut", bbox_inches='tight', dpi=600)

    err_e = np.log10(np.absolute((enr[:,1] - enr[0,1]) / enr[0,1]))

    ax[1].plot(enr[:,0], err_e)
    ax[1].tick_params(which="both", direction="in")
    ax[1].set_xlabel("time [kyr]")
    ax[1].set_ylabel(r"$\Delta_{10} E$")
    pt.show()
    
    
if __name__ == "__main__":
    PlotOutput()
    
