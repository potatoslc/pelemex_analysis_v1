"""
Author: u0890475
potatoslc version of PV_gc.py
Need to rewrtie most of the code


"""
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pandas as pd
import re
import os
from all_functions import *
"""
# define the file path over here
"""
#path = "/media/u1474458/Data/phi375_amr3_plt/plt/" previous file path
#path ="/media/u0890475/6f7f5b18-6951-4d23-9e1b-146d3d4c2671/Matt_Old/phi375_amr3_plt/plt/"
path ='/media/u0890475/6f7f5b18-6951-4d23-9e1b-146d3d4c2671/Matt_Old/phi375_amr3_plt/plt/'
U_title = "U_wall_phi375" # name of all CSV files before timestep
Flame_title = "CH2_0.8_phi375" # name of all CSV files before timestep


"""
change if needed
For the current case, the names of x,y coord are:
x-coordinate: Points:0
y-coordinate: Points:1
"""


"""
get all files with needed prefix
also filter out file with errors

combine file names with same series number
as:[0, 'U_wall_phi375_0.csv', 'CH2_0.8_phi375_0.csv']
"""

speed_files = file_findall_complex(U_title, path)
flame_files = file_findall_complex(Flame_title, path)

speed_dict = dict(speed_files)
flame_dict = dict(flame_files)
match_both = [ [x,speed_dict[x],flame_dict[x]] for x in speed_dict.keys() ]
match_both.sort(key=lambda x:x[0])

num_ttl_file = len(match_both)

result_array = []
for i in range(num_ttl_file):
    paradata_filename = path+match_both[i][1]
    isosurface_filename = path+match_both[i][2]
    result_data = get_v_gc(paradata_filename, isosurface_filename)
    print("finishing process data: "+str(i)+ ". And g is: "+str(result_data[0]))
    result_array.append(result_data)

result_df = pd.DataFrame(result_array,columns = ["g_isosrf","x_minx_isosrf","y_minx_isosrf","x_minx_para","y_minx_para","Timestep"])
result_df.to_csv("gc_time.csv",index=False)




"""
dc = Sl_ext/g (dc should be an array) sl:from cantera,g:velocity gradient
"""


"""
dt = sqrt(thermal_diffusivity/g) thermal_diffusivity:cantera, g:velocity gradient
"""
    



