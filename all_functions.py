#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 23 11:28:59 2025

@author: u0890475
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import h5py
import time
import cantera as ct

import re
def derivative(x, y):
    if len(x)!=len(y):
        print("Input files have different size")
        return
    

    dx = np.diff(x)
    dy = np.diff(y)
    dydx = [y/x if x!= 0 else 0 for (y,x) in zip(dy,dx) ]
    #dydx[0:-1] = dy/dx
    
    dydx.append((y[-1] - y[-2])/(x[-1] - x[-2]))

    return dydx
def file_findall_complex(header,directory):
    file_seq = []
    file_match_patten = header+r"_(\d+)\.csv"
    file_list = os.listdir(directory)
    uin_files =[x for x in file_list if header in x]
    total_files = len(uin_files)

    
    for i in range(total_files):
        temp_filename = uin_files[i]
        match1= re.search(file_match_patten,temp_filename)
        if(match1):
            file_seq.append([int(match1.group(1)),temp_filename])
            
        else:
            continue
    
    file_seq.sort(key=lambda x:x[0])
    return file_seq


def get_v_gc(paradata_filename,isosurface_filename):
    paradata_df = pd.read_csv(paradata_filename)
    isosurface_df = pd.read_csv(isosurface_filename)
    current_time = paradata_df["Time"][0]
    """
    get u info from isosurface file
    x: Points:0
    y: Points:1
    x velocity: x_velocity
    """
    minx_loc = isosurface_df["Points:0"].argmin()
    x_minx ,y_minx, xvel_minx= isosurface_df.loc[minx_loc,["Points:0"]],isosurface_df.loc[minx_loc,["Points:1"]],isosurface_df.loc[minx_loc,["x_velocity"]]




    """
    get speed profile from paraview file
    """
    flame_minx_loc = paradata_df["Points:0"].argmin()
    x_flame_minx,y_flame_minx,z_flame_minx= paradata_df.loc[flame_minx_loc,["Points:0"]],paradata_df.loc[flame_minx_loc,["Points:1"]],paradata_df.loc[flame_minx_loc,["Points:2"]]
    #x_vel_flame_minx=paradata_df.loc[flame_minx_loc,["x_velocity"]]
    x_vel_filter = np.where(paradata_df["x_velocity"]<=0,1,0)
    paradata_df["x_filter"] = x_vel_filter*paradata_df["Points:0"]
    #paradata_df["y_filter"] = x_vel_filter*paradata_df["Points:1"]
    #paradata_df["z_filter"] = x_vel_filter*paradata_df["Points:2"]
    x_fb = min(paradata_df[paradata_df["x_filter"]>0]["Points:0"])
    
    len_sep = max(paradata_df["x_filter"] - x_fb)
    #x_diff = abs(paradata_df["Points:0"] - (1.25*len_sep - x_minx))
    x_comp = (1.25*len_sep - x_minx)
    #t1= time.perf_counter()
    x_diff = [ (x-x_comp ) for x in paradata_df["Points:0"]]
    #t2 = time.perf_counter()
    #print("takes:" +str(t2-t1))
    min_diff_loc = np.array(x_diff).argmin()
    dv = float(paradata_df["x_velocity"][min_diff_loc])
    dz = float(paradata_df["Points:1"][min_diff_loc])
    g = dv/dz

    return [g, float(x_flame_minx),float(y_flame_minx),float(x_minx),float(y_minx),current_time]


def build_flame(T,P,phi,fuel='H2:1',vel=4):
    loglevel = 1
    g = ct.Solution('mechanism.yaml')
    mdotR = vel*g.density
    ct.transport_model = 'multicomponent'
    air = "O2:0.209,N2:0.791" # mol fraction of air with only N2,O2
    g.set_equivalence_ratio(phi=phi, fuel="H2:1", oxidizer=air)
	#g.X = {'H2':1, 'O2':.5/phi, 'N2':(.5*3.76)/phi}
    g.TP = T,P
    gridd = np.linspace(0,.01,1024)
    flame = ct.FreeFlame(gas = g,grid = gridd)
    #https://cantera.org/2.0/doxygen/html/classCantera_1_1FreeFlame.html
    flame.transport_model = 'multicomponent'
    flame.inlet.mdot = mdotR
    flame.set_refine_criteria(ratio=4, slope=0.05, curve=0.1, prune=0.01)
    flame.solve(loglevel=loglevel, refine_grid = False, auto = False)

    return flame


def get_flame_info(flame):
    try:
        flame_grid = flame.grid
        speed = flame.velocity
        
        
    except:
        print("flame has erros")
        return 

    strain_rate = derivative(flame_grid, speed)
    max_loc = np.argmax(strain_rate)
    thermal_diffus = (flame.thermal_conductivity[max_loc])/(flame.cp_mass[max_loc])/(flame.density[max_loc])
    sd = flame.velocity[-1]

    return sd,thermal_diffus



def get_gc_data(path):
    try:
        df = pd.read_csv(path)
    except:
        return

    return df

def mono_plot(x,y,xlabel="",ylabel="",title=""):
    fig,ax = plt.subplots()
    ax.plot(x,y,'r')
    ax.set(xlabel=xlabel,ylabel=ylabel,title=title)

    return plt

def duo_plot_one_scale(x1,y1,x2,y2,label_list= [" "," "," "]):
    labels_ct = len(label_list)
    x1_label,y1_label,y2_label =label_list[0],label_list[1],label_list[2]
    if labels_ct ==2:
        x2_label = x1_label
        y2_label = y1_label
    if labels_ct==3:
        y2_label = label_list[2]    
    plt.plot(x1,y1 ,label=y1_label) #,label=y1_label
    plt.plot(x2,y2 ,label=y2_label) # ,label=y2_label
    plt.legend()
    plt.xlabel(x1_label)
    return plt

def duo_plot_duo_scale(x1,y1,x2,y2,label_list= [" "," "," "]):
    x1_label,y1_label,y2_label =label_list[0],label_list[1],label_list[2]
    labels_ct = len(label_list)
    if labels_ct ==2:
        x2_label = x1_label
        y2_label = y1_label
    if labels_ct==3:
        y2_label = label_list[2]

    fig, ax1 = plt.subplots()
    print(y1_label,y2_label)
    ax1.set_xlabel(x1_label)
    ax1.set_ylabel(y1_label)
    ax1.legend([y1_label,y2_label])
    ln1 = ax1.plot(x1,y1,'r',label=y1_label)
    #ax1.tick_params(axis='y',labelcolor='r')

    ax2 = ax1.twinx()
    ax2.set_ylabel(y2_label)
    ln2 = ax2.plot(x2,y2,'b',label = y2_label)
    leg = ln1 + ln2
    labels = [x.get_label() for x in leg]
    ax1.legend(leg,labels,loc=0)
    #ax2.legend([y2_label])
    #ax2.legend()
    #ax2.tick_params(axis='y', labelcolor='b')
    #fig.legend()
    #plt.legend()


    return plt

