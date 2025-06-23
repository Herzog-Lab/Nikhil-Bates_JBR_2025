# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 11:18:48 2021

@author: Nikhil
"""

exec(open('./Imports/modules_pixar.py').read())

path = input("Enter full name of the results folder:")
cache = pd.read_csv(str(path) + "/cache/cache.csv", index_col = 0)
#creates folders to save results
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

createFolder('./' +str(path) + '/Replotted Data')
createFolder('./' +str(path) + '/Replotted Data')
createFolder('./' +str(path) + '/Replotted Data/Maps')
createFolder('./' +str(path) + '/Replotted Data/Maps/Amp Maps')
createFolder('./' +str(path) + '/Replotted Data/Maps/Period Maps')
createFolder('./' +str(path) + '/Replotted Data/Maps/Phase Maps')



sampling_interval_mins = cache.iloc[0,0]
min_per = cache.iloc[0,1]
max_per = cache.iloc[0,2]

if phasemap_time == 'auto':
    phasemap_time = cache.iloc[0,3]
else:
    phasemap_time = phasemap_time
compiled_phase = pd.read_csv(str(path) + "/Compiled Phase.csv")
compiled_period = pd.read_csv(str(path) + "/Compiled Period.csv")
cwt_pha_data = pd.read_csv(str(path) + "/CWT Instantaneous Phase_rad.csv", index_col = 0)
cwt_pha_data = cwt_pha_data.T.reset_index(drop = True)
#cwt_pha_data.columns = np.arange(0 , len(cwt_pha_data.columns)*sampling_interval_mins,sampling_interval_mins)

heatmap_data = pd.read_csv(str(path) + "/Smoothed Traces.csv", index_col = 0)


time = np.arange(0 , len(cwt_pha_data.columns)*1.0*sampling_interval_mins,sampling_interval_mins)                        
while phasemap_time*60  not in time:
    warnings.warn("Time entered for phasemap generation is either out of range or not compatible with sampling interval. Please enter a different time")
    phasemap_time = float(input("Enter time (IN HOURS) for phasemap generation : "))


#calculates relative phases and amp and phase at time t
phaset_meta = ((compiled_phase["Metacyc Pha (t = 0h)"]*2*np.pi/compiled_period["Metacyc Per (h)"]) + 
               (phasemap_time* 2 * np.pi/compiled_period["Metacyc Per (h)"]))%(2*np.pi)
               
meta_meanphase = (np.arctan2(np.sin(phaset_meta).sum(),np.cos(phaset_meta).sum()))%(2*np.pi)
meta_rel_phase = phaset_meta - meta_meanphase

for i in np.arange(0,len(meta_rel_phase)):
    if meta_rel_phase[i] > np.pi:
        meta_rel_phase[i] = meta_rel_phase[i]%-np.pi
    elif meta_rel_phase[i] < -1* np.pi:
        meta_rel_phase[i] =  meta_rel_phase[i]%np.pi
    else: continue

t_col = int(phasemap_time* 60/sampling_interval_mins)
phaset_cwt = cwt_pha_data[t_col]

cwt_meanphase = (np.arctan2(np.sin(phaset_cwt).sum(),np.cos(phaset_cwt).sum()))%(2*np.pi)
cwt_rel_phase =phaset_cwt - cwt_meanphase
for i in np.arange(0,len(cwt_rel_phase)):
    if cwt_rel_phase[i] > np.pi:
        cwt_rel_phase[i] = cwt_rel_phase[i]%-np.pi
    elif cwt_rel_phase[i] < -1* np.pi:
        cwt_rel_phase[i] = cwt_rel_phase[i]%np.pi
    else: continue



compiled_amplitude = pd.read_csv(str(path) + "/Compiled Amplitude.csv") 
cwt_meanamp = round(compiled_amplitude['Avg CWT Amp'].mean(),5)
meta_meanamp = round(compiled_amplitude['Metacyc Amp'].mean(),5)
xov_meanamp =  round(compiled_amplitude['XOV Amp'].mean(),5)
meta_meanphaseh = round(meta_meanphase*12/np.pi,)
cwt_meanphaseh = round(cwt_meanphase*12/np.pi,1)
meta_meanper = round(compiled_period["Metacyc Per (h)"].mean(),2)
cwt_meanper = round(compiled_period["Avg CWT Per (h)"].mean(), 2)
FFT_meanper = round(compiled_period["FFT Per (h)"].mean(),2)
pslope_meanper = round(compiled_period['XOV + Slope Per (h)'].mean(),2)
nslope_meanper = round(compiled_period['XOV - Slope Per (h)'].mean(),2)
peak_meanper = round(compiled_period['XOV Peak Per (h)'].mean(),2)
trough_meanper = round(compiled_period['XOV Trough Per (h)'].mean(),2)


#sets min and max amplitudes
if cwt_minamp == 'auto' and cwt_maxamp == 'auto':
    cwt_minamp = compiled_amplitude['Avg CWT Amp'].min()
    cwt_maxamp = compiled_amplitude['Avg CWT Amp'].max()


if metacyc_minamp == 'auto' and metacyc_maxamp == 'auto':
    metacyc_minamp = compiled_amplitude['Metacyc Amp'].min()
    metacyc_maxamp = compiled_amplitude['Metacyc Amp'].max()


if xov_minamp == 'auto' and xov_maxamp == 'auto':
    xov_minamp = compiled_amplitude['XOV Amp'].min()
    xov_maxamp = compiled_amplitude['XOV Amp'].max()


#sets min and max pers

xov_period = pd.concat([compiled_period['XOV + Slope Per (h)'], compiled_period['XOV - Slope Per (h)'],
                        compiled_period['XOV Peak Per (h)'], compiled_period['XOV Trough Per (h)']], axis = 1 )




def mapper(z,vmin,vmax,colormap,title,cbarlabel, bins, tick):
    x = compiled_period['X']
    y = compiled_period['Y']
    y_vals, y_idx = np.unique(y, return_inverse=True)
    x_vals, x_idx = np.unique(x, return_inverse=True)
    vals_array = np.empty(y_vals.shape + x_vals.shape)
    vals_array.fill(np.nan) # or whatever yor desired missing data flag is
    vals_array[y_idx, x_idx] = z
    #extent = x_min, x_max, y_min, y_max = [min(x)-10, max(x)+ 10, min(y)-10,max(y)+10]
    fig, ax = plt.subplots(figsize=(4,4.5), tight_layout = True)
    ax.set_aspect(1)
    #ax.set_xlim(min(x)-2, max(x)+2)
    #ax.set_ylim(min(y)-2, max(y)+2)
    ax.invert_yaxis()
    #ax.invert_xaxis()
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    ax.set_title(title, fontsize ="13", pad = 8)
    ax.grid(c='k', alpha=0)
    cm = plt.cm.get_cmap(colormap)
    normalize = mpl.colors.Normalize(vmin= vmin, vmax= vmax)
    sc = plt.imshow(vals_array,cmap=cm, norm = normalize)
    ax.set_facecolor(map_background_color)
    if tick == "amplitude":
        ticks = np.linspace(vmin, vmax, bins)
    elif tick =="period":
        ticks = np.linspace(vmin, vmax,7)
    elif tick == "phase":
        ticks = np.arange(vmin, vmax+1, 3)
    cbar = plt.colorbar(sc, fraction=0.046, pad=0.04, ticks = ticks)
    cbar.set_label(cbarlabel, size = '11', rotation=90)
    cbar.ax.tick_params(labelsize=10)

if maps_minper == 'auto' and maps_maxper  == 'auto':
    vmin_meta = round(compiled_period["Metacyc Per (h)"].min(),1)
    vmax_meta = round(compiled_period["Metacyc Per (h)"].max() ,1)
    vmin_cwt =round(compiled_period["Avg CWT Per (h)"].min(),1)
    vmax_cwt = round(compiled_period["Avg CWT Per (h)"].max(),1)     
    vmin_xov = round(xov_period.min().min(),1)
    vmax_xov = round(xov_period.max().max(),1)    
    vmin_fft = round(compiled_period["FFT Per (h)"].min(),1)
    vmax_fft = round(compiled_period['FFT Per (h)'].max(),1) 
else:
    vmin_meta = vmin_cwt = vmin_fft = vmin_xov = maps_minper 
    vmax_meta = vmax_cwt = vmax_fft = vmax_xov = maps_maxper 
 
    
#CWT amplitude map
mapper(compiled_amplitude['Avg CWT Amp'],cwt_minamp,cwt_maxamp,
       periodmap_color,'Avg CWT Amp ' r'$\mu$' + ': ' + str(cwt_meanamp),'amplitude', 8, tick = "amplitude")
plt.savefig(str(path) + "/Replotted Data/Maps/Amp Maps/Avg CWT Amp Map .tif", format="TIF", dpi=dpi)    
 
#metacycle amplitude map
mapper(compiled_amplitude['Metacyc Amp'],metacyc_minamp,metacyc_maxamp,
       periodmap_color,'Metacyc Amp ' r'$\mu$' +': ' + str(meta_meanamp),'amplitude', 8, tick = "amplitude")
plt.savefig(str(path) + "/Replotted Data/Maps/Amp Maps/Metacyc rAmp Map .tif", format="TIF", dpi=dpi)

#XOV amplitude map
mapper(compiled_amplitude['XOV Amp'],xov_minamp, xov_maxamp,
       periodmap_color,'Crossover Amp ' r'$\mu$' + ': ' + str(xov_meanamp),'amplitude', 8, tick = "amplitude")
plt.savefig(str(path) + "/Replotted Data/Maps/Amp Maps/Crossover Amp Map .tif", format="TIF", dpi=dpi)   


#metacycle period map
mapper(compiled_period["Metacyc Per (h)"],vmin_meta,vmax_meta, periodmap_color, 
       'Metacyc Per (' r'$\mu$' +': ' +str(meta_meanper) +'h)','period (h)', 8, tick = "period")
plt.savefig(str(path) + "/Replotted Data/Maps/Period Maps/Metacyc Period Map .tif", format="TIF", dpi=dpi)

#CWT period map
mapper(compiled_period["Avg CWT Per (h)"],vmin_cwt,vmax_cwt, periodmap_color, 'CWT Per ('r'$\mu$' +': ' +str(cwt_meanper) + 'h)',
       'period (h)', 8,tick = "period")
plt.savefig(str(path) + "/Replotted Data/Maps/Period Maps/CWT Period Map .tif", format="TIF", dpi=dpi)

#Fourier period map
mapper(compiled_period['FFT Per (h)'],vmin_fft,vmax_fft, periodmap_color, 'FFT Per ('r'$\mu$' +': ' +str(FFT_meanper) +'h)',
 'period (h)', 8, tick = "period")
plt.savefig(str(path) + "/Replotted Data//Maps/Period Maps/FFT Period Map .tif", format="TIF", dpi=dpi)

#posslope period map
mapper(compiled_period["XOV + Slope Per (h)"],vmin_xov, vmax_xov, periodmap_color, 
       '+ Slope Per ('r'$\mu$' +': ' +str(pslope_meanper) + 'h)','period (h)', 8, tick = "period")
plt.savefig(str(path) + "/Replotted Data/Maps/Period Maps/Pos Slope Period Map .tif", format="TIF", dpi=dpi)

#negslope period map
mapper(compiled_period["XOV - Slope Per (h)"],vmin_xov, vmax_xov, periodmap_color, 
       '- Slope Per ('r'$\mu$' +': ' +str(nslope_meanper) + 'h)','period (h)', 8, tick = "period")
plt.savefig(str(path) + "/Replotted Data/Maps/Period Maps/Neg Slope Period Map .tif", format="TIF", dpi=dpi)

#peak period map
mapper(compiled_period["XOV Peak Per (h)"],vmin_xov, vmax_xov, periodmap_color, 
       'Peak Per ('r'$\mu$' +': ' +str(peak_meanper) + 'h)','period (h)', 8, tick = "period")
plt.savefig(str(path) + "/Replotted Data/Maps/Period Maps/Peak Period Map .tif", format="TIF", dpi=dpi)

#trough period map
mapper(compiled_period["XOV Trough Per (h)"],vmin_xov, vmax_xov, periodmap_color, 
       'Trough Per ('r'$\mu$' +': ' +str(trough_meanper) + 'h)','period (h)', 8, tick = "period")
plt.savefig(str(path) + "/Replotted Data/Maps/Period Maps/Trough Period Map .tif", format="TIF", dpi=dpi)


#metacycle rel phase map
mapper(meta_rel_phase *12/np.pi,-12,12, 
       phasemap_color, 'rPha at ' + str(phasemap_time) + ' h (Metacyc) 'r'$\mu$' +': '  
       + str(meta_meanphaseh) + 'h', 'phase diff (h)', 8, tick = "phase")
plt.savefig(str(path) + "/Plots/Maps/Phase Maps/Metacyc rPhase Map .tif", format="TIF", dpi=dpi)

#metacycle abs phase map
mapper(phaset_meta*12/np.pi,0,24, 
       phasemap_color, 'Pha at ' + str(phasemap_time) + 'h (Metacyc) 'r'$\mu$' +': '  
       + str(meta_meanphaseh) + 'h','phase (h)', 8, tick = "phase")
plt.savefig(str(path) + "/Replotted Data/Maps/Phase Maps/Metacyc Phase Map .tif", format="TIF", dpi=dpi)

#CWT Rel phase map
mapper(cwt_rel_phase*12/np.pi,-12,12, 
       phasemap_color, 'rPhase at ' +str(phasemap_time) + 'h (CWT) 'r'$\mu$' +': '  
       + str(cwt_meanphaseh) +'h', 'phase diff (h)', 8, tick = "phase")

plt.savefig(str(path) + "/Replotted Data/Maps/Phase Maps/CWT rPhase Map .tif", format="TIF", dpi=dpi)

#CWT Abs phase map
mapper(phaset_cwt *12/np.pi, 0,24, 
       phasemap_color, 'Phase at ' +str(phasemap_time) + ' h (CWT) ' r'$\mu$' +': '  
       + str(cwt_meanphaseh) +'h', 'phase (h)', 8, tick = "phase")

plt.savefig(str(path) + "/Replotted Data/Maps/Phase Maps/CWT Phase Map.tif", format="TIF", dpi=dpi)




compiled_phase.iloc[:,4] = round(phaset_meta*12/np.pi, 2)
compiled_phase.iloc[:,5] = round(meta_rel_phase*12/np.pi, 2)
compiled_phase.iloc[:,7] = round(phaset_cwt*12/np.pi,2)
compiled_phase.iloc[:,8] = round(cwt_rel_phase*12/np.pi,2)
compiled_phase.columns = ["X","Y","CycID/Pixel","Metacyc Pha (t = 0h)","Metacyc Pha (t = " + str(phasemap_time)+"h)",	
                          "Metacyc rPha (t = " + str(phasemap_time)+"h)",
                          "CWT Pha (t = 0h)","CWT Pha (t = " + str(phasemap_time)+"h)","CWT rPha (t = " + str(phasemap_time)+"h)",
                          'XOV + Slope Pha (h)',
                           "XOV + Slope rPha (h)","XOV - Slope Pha (h)","XOV - Slope rPha (h)",	"XOV Peak Pha (h)",
                           "XOV Peak rPha (h)","XOV Trough Pha (h)","XOV Trough rPha (h)"]

compiled_phase.to_csv("./" +str(path) + "/Replotted Data/Replotted Phase Data.csv", index = False)

#Plots heatmap of CWT phases. 'method' and 'metric' can be changed based on user requirements for better plots
#https://docs.scipy.org/doc/scipy/reference/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
#https://docs.scipy.org/doc/scipy/reference/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
sns.clustermap(heatmap_data, method=cluster_method, metric=distance_metric,
               z_score = z_score, row_cluster=True, col_cluster=False, row_linkage=None, standard_scale = None,
               cmap = clustermap_color, cbar_pos = None)
plt.savefig("./" +str(path) + "/Replotted Data/HeatMap.tif", format = 'TIF', dpi=dpi)



print("Replotting Completed. Results saved in Replotted Data folder inside the Results folder")