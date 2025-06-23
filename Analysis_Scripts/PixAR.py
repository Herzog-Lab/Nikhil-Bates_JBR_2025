# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
#------------------------------------------------------------------------------
# Analysis Parameters (modify if required)
min_per = 18  #min period (in h) to be analysed
max_per = 30  #max period (in h) to be analysed
pixels_per_image = 1000 #use lower value to reduce analysis time (if required). 
                        #Note: depending on your data dpi, lower pixels_per_image values may lead to reduced spatial resolution 
file_extention = ".tif"
nper = 200  #number of periods to analyse between min and max per (for CWT analysis only)
T_c = 48  #cutoff period (in h) for sinc filter. 48 works for most cases and 
            #may need to be changed ONLY under some circumstances.
smooth = 'default' #Extent of smoothing (sigma for gaussian blur) for metacycle and cwt analysis. 
                #Default is 1 and works for most cases. But if the traces look extremely noisy even after this, 
                #change to values to 2 or 3 (without quotes). Unlikely that you will have to modify this or go higher than 3. 
   
                
#------------------------------------------------------------------------------
#Plot Customization (modify if required)
#List of colormaps for maps: https://matplotlib.org/stable/tutorials/colors/colormaps.html
#List of colors for plots: https://matplotlib.org/stable/gallery/color/named_colors.html
#set colors
periodmap_color = "jet" #colorscale for periodmaps
phasemap_color = "hsv" #colorscale for phasemaps (preferably use cyclic colormaps)
map_background_color = 'black'#background color for all plots mentioned above
clustermap_color = "summer_r" #colormaps for clustermap

#set min and max periods for period maps
maps_minper = 'auto'
maps_maxper = 'auto'


#set dpi for saving figures
dpi = 300
#-------------------------------------------------------------------------------
#Do not change anything in here
exec(open('./Imports/modules_pixar.py').read())
exec(open('./Imports/scaling_factor.py').read())
exec(open('./Imports/metacyc_pixar.py').read())
print("MetaCycle Analysis Completed. CWT Analysis in progress...")
exec(open('./Imports/cwt_pixar.py').read())
print("CWT Analysis Completed. Crossover Analysis in progress...")
exec(open('./Imports/crossover_pixar.py').read())
print("Crossover Analysis Completed. Saving Data and Generating Plots...")
exec(open('./Imports/plotter_pixar.py').read())
