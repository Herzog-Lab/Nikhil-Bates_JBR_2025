#saves arrhythmic time series
ar_pix_val_ts.to_csv("./Results_" +str(filename_raw)+ "/Arrhythmic Traces.csv")

#saves smoothed Traces
pd.DataFrame(xov_smoothed).to_csv("./Results_" +str(filename_raw)+ "/Smoothed Traces.csv")

#plots  randomly selected time series used for crossover analysis (sigma = 3)
if len(xov_smoothed) > 100:
    n =100
else:
    n = int(len(xov_smoothed)/2)

for i in range(n):
    sample = xov_smoothed[random.randint(0, len(xov_smoothed)-1)]
    plt.plot(xticks[0], sample)
plt.xlabel ("time (h)")
plt.ylabel("Normalized Pixel Intensity") 
plt.title(" Crossover Analysis Representative Time Series (sigma =3)")  
for k in vlines:
    plt.axvline(x=k, color='k', linestyle='--', alpha = 0.2)

plt.savefig("./Results_" +str(filename_raw)+ "/Plots/Sample timeseries(smoothed).tif", format="TIF", dpi=dpi) 
plt.show()
plt.close()


#plots  phases  for randomly selected pixel
for i in range(n):
    sample = cwt_phase.loc[random.randint(0, len(cwt_phase)-1)]
    plt.plot(sample*12/np.pi)
plt.title("Randomly Sampled CWT Phase traces")
plt.xlabel("Time (h)")
plt.ylabel("Phase (h)")
plt.savefig('./Results_' +str(filename_raw)+ '/Plots/Sample CWT Phases.tif', format="TIF", dpi=dpi)
plt.show()
plt.close()

#saves crossover data
xov_amp_cycwise.to_csv("./Results_" +str(filename_raw)+"/Crossover cyclewise amp.csv")
#xov_cycwise_per.to_csv("./Results_" +str(filename_raw)+"/Crossover cycwise per.csv")
xov_posslope_period.to_csv("./Results_" +str(filename_raw)+"/Crossover pos slope per.csv")
xov_negslope_period.to_csv("./Results_" +str(filename_raw)+"/Crossover neg slope per.csv")
peaks_period.to_csv("./Results_" +str(filename_raw)+"/Crossover peak per.csv")
troughs_period.to_csv("./Results_" +str(filename_raw)+ "/Crossover trough per.csv")
xov_period = pd.concat([xov_posslope_period.mean(axis=1),xov_negslope_period.mean(axis = 1),
                          peaks_period.mean(axis = 1),troughs_period.mean(axis = 1)], axis = 1)
xov_period.columns = ['XOV + Slope Per (h)', 'XOV - Slope Per (h)', 'XOV Peak Per (h)', 'XOV Trough Per (h)']
xov_period.to_csv("./Results_" +str(filename_raw)+ "/Crossover Periods.csv")

cycwise_posslope_phase.to_csv("./Results_" +str(filename_raw)+ "/Crossover pos slope phase.csv")
cycwise_negslope_phase.to_csv("./Results_" +str(filename_raw)+ "/Crossover neg slope phase.csv")
cycwise_peak_phase.to_csv("./Results_" +str(filename_raw)+ "/Crossover peak phase.csv")
cycwise_trough_phase.to_csv("./Results_" +str(filename_raw)+ "/Crossover trough phase.csv")

#saves compiled analysis results
avg_cwt_per = avg_cwt_per.reset_index(drop = True)
compiled_period = pd.concat([rhythmic_pixs['X'],rhythmic_pixs['Y'],rhythmic_pixs['CycID'], rhythmic_pixs['meta2d_BH.Q'], 
                       rhythmic_pixs['meta2d_period'], round(avg_cwt_per,2), round(FFT_period,2), xov_period], axis = 1)
compiled_period.columns=['X', 'Y', 'CycID/Pixel','Metacyc BH.Q', 'Metacyc Per (h)','Avg CWT Per (h)', 'FFT Per (h)',
                         'XOV + Slope Per (h)', 'XOV - Slope Per (h)', 'XOV Peak Per (h)', 'XOV Trough Per (h)']



compiled_amplitude = pd.concat([rhythmic_pixs['X'],rhythmic_pixs['Y'],rhythmic_pixs['CycID'], rhythmic_pixs['meta2d_rAMP'],avg_cwt_amp,
                                xov_amp_cycwise.mean(axis = 1)], axis = 1) 
compiled_amplitude.columns=['X', 'Y', 'CycID/Pixel', 'Metacyc Amp','Avg CWT Amp', 'XOV Amp']

                                
compiled_phase = pd.concat([rhythmic_pixs['X'],rhythmic_pixs['Y'],rhythmic_pixs['CycID'], rhythmic_pixs['meta2d_phase']], axis = 1)
compiled_phase.columns = [ 'X', 'Y', 'CycID/Pixel', 'Metacyc Pha (t = 0h)']



#plots com on top of tissue at user defined time
time_range = pd.DataFrame(time).T
t = int(phasemap_time * (60/sampling_interval_mins))
if image[t].max() > 300:
    normalize = mpl.colors.Normalize(vmin= 300, vmax= image[t].max())
else:
    normalize = mpl.colors.Normalize(vmin= image[t].min(), vmax= image[t].max())
plt.imshow(image[t], cmap = 'viridis', norm = normalize )
plt.scatter(com.loc[t][1],com.loc[t][0],color = 'red', s = 10)
plt.axis('off')
plt.title(label = "COM (t = "+ str(phasemap_time) + "h)", fontsize = "10")
plt.savefig("./Results_" +str(filename_raw)+ "/Plots/COM.tif", format="TIF", dpi=dpi)

for i in range(len(time_range)):
    #print(i)
    #normalize = mpl.colors.Normalize(vmin=300, vmax= image[i].max())
    plt.imshow(image[i], cmap = 'viridis', norm = normalize)
    plt.scatter(com.loc[i][1],com.loc[i][0],color = 'red', s = 10)
    plt.axis('off')
    plt.title(label = "COM (t = "+ str(round(time_range[0][i]/60,2)) + "h)", fontsize = "10")
    plt.savefig("./Results_" +str(filename_raw)+"/Plots/Stack Plots/COM/COM(t = " + str(round(time_range[0][i]/60,2)) + "h).tif", format="TIF", dpi=dpi)
    plt.close('all')


#plots COM 
ncycles_com = int(len(com)*sampling_interval_mins/(24*60))
com_time = np.arange(0, int((24*ncycles_com*(sampling_interval_mins/60))+1), int(24*sampling_interval_mins/60))
fig, ax = plt.subplots(figsize=(4,4.5), tight_layout = True)
legend =[]
for i in np.arange(1, ncycles_com+2):
    legend.append('day ' + str(i))
for i in com_time:
    plt.plot(com['com_y'].loc[i : i + int(24*sampling_interval_mins/60)], 
             com['com_x'].loc[i: i+int(24*sampling_interval_mins/60)], linestyle='dashed', marker='o', 
             linewidth=1.25, markersize = 3.25)
    ax.invert_yaxis()
    ax.set_ylabel("Y Pixels", fontsize = "15")
    ax.set_xlabel("X Pixels", fontsize = "15")
    #ax.set_yticks([])
    #ax.set_xticks([])
    ax.legend(legend, prop={'size': 5})
    ax.set_title("COM Trajectory", fontsize = "18")
    
plt.savefig("./Results_" +str(filename_raw)+"/Plots/COM_Trajectory.tif", format="TIF", dpi=dpi)
plt.show()
plt.close()




#calculates relative phases at time t
phaset_meta = ((compiled_phase["Metacyc Pha (t = 0h)"]*2*np.pi/compiled_period["Metacyc Per (h)"]) 
               + (phasemap_time* 2 * np.pi/compiled_period["Metacyc Per (h)"]))%(2*np.pi)
               
meta_meanphase = (np.arctan2(np.sin(phaset_meta).sum(),np.cos(phaset_meta).sum()))%(2*np.pi)
meta_rel_phase = phaset_meta - meta_meanphase

for i in range(len(meta_rel_phase)):
    if meta_rel_phase[i] > np.pi:
        meta_rel_phase[i] = meta_rel_phase[i]%-np.pi
    elif meta_rel_phase[i] < -1* np.pi:
        meta_rel_phase[i] =  meta_rel_phase[i]%np.pi
    else: continue

phaset_cwt = pd.DataFrame(cwt_phase)[phasemap_time* 60/sampling_interval_mins]
cwt_meanphase = (np.arctan2(np.sin(phaset_cwt).sum(),np.cos(phaset_cwt).sum()))%(2*np.pi)
cwt_rel_phase = phaset_cwt - cwt_meanphase
for i in range(len(cwt_rel_phase)):
    if cwt_rel_phase[i] > np.pi:
        cwt_rel_phase[i] = cwt_rel_phase[i]%-np.pi
    elif cwt_rel_phase[i] < -1* np.pi:
        cwt_rel_phase[i] = cwt_rel_phase[i]%np.pi
    else: continue

cwt_meanamp = round(compiled_amplitude['Avg CWT Amp'].mean(),5)
meta_meanamp = round(compiled_amplitude['Metacyc Amp'].mean(),5)
xov_meanamp =  round(compiled_amplitude['XOV Amp'].mean(),5)
meta_meanphaseh = round(meta_meanphase*12/np.pi,2)
cwt_meanphaseh = round(cwt_meanphase*12/np.pi,2)
meta_meanper = round(compiled_period["Metacyc Per (h)"].mean(),2)
cwt_meanper = round(compiled_period["Avg CWT Per (h)"].mean(), 2)
FFT_meanper = round(compiled_period["FFT Per (h)"].mean(),2)
pslope_meanper = round(compiled_period['XOV + Slope Per (h)'].mean(),2)
nslope_meanper = round(compiled_period['XOV - Slope Per (h)'].mean(),2)
peak_meanper = round(compiled_period['XOV Peak Per (h)'].mean(),2)
trough_meanper = round(compiled_period['XOV Trough Per (h)'].mean(),2)

compiled_phase.insert(4, "Metacyc Pha (t = " + str(phasemap_time)+"h)", round(phaset_meta*12/np.pi, 2))
compiled_phase.insert(5, "Metacyc rPha (t = " + str(phasemap_time)+"h)", round(meta_rel_phase*12/np.pi, 2))
compiled_phase.insert(6, "CWT Pha (t = 0h)", round(cwt_phase[0]*12/np.pi,2))
compiled_phase.insert(7, "CWT Pha (t = "+ str(phasemap_time)+"h)", round(phaset_cwt*12/np.pi,2))
compiled_phase.insert(8, "CWT rPha (t = " + str(phasemap_time)+"h)", round(cwt_rel_phase*12/np.pi,2))
compiled_phase = pd.concat([compiled_phase, xov_phase], axis = 1)
compiled_amplitude.to_csv("./Results_" +str(filename_raw)+"/Compiled Amplitude.csv", index = False)
compiled_period.to_csv("./Results_" +str(filename_raw)+"/Compiled Period.csv", index = False)
compiled_phase.to_csv("./Results_" +str(filename_raw)+"/Compiled Phase.csv", index = False)

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
    ax.invert_xaxis()
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
mapper(compiled_amplitude['Avg CWT Amp'],compiled_amplitude['Avg CWT Amp'].min(),compiled_amplitude['Avg CWT Amp'].max(),
       periodmap_color,'Avg CWT Amp ' r'$\mu$' + ': ' + str(cwt_meanamp),'amplitude', 8, tick = "amplitude")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Amp Maps/Avg CWT Amp Map .tif", format="TIF", dpi=dpi)    
 
#metacycle amplitude map
mapper(compiled_amplitude['Metacyc Amp'],compiled_amplitude['Metacyc Amp'].min(),compiled_amplitude['Metacyc Amp'].max(),
       periodmap_color,'Metacyc Amp ' r'$\mu$' +': ' +str(meta_meanamp),'amplitude', 8, tick = "amplitude")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Amp Maps/Metacyc Amp Map .tif", format="TIF", dpi=dpi)

#XOV amplitude map
mapper(compiled_amplitude['XOV Amp'],compiled_amplitude['XOV Amp'].min(),compiled_amplitude['XOV Amp'].max(),
       periodmap_color,'Crossover Amp ' r'$\mu$' + ': ' + str(xov_meanamp),'amplitude', 8, tick = "amplitude")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Amp Maps/Crossover Amp Map .tif", format="TIF", dpi=dpi)   


#metacycle period map
mapper(compiled_period["Metacyc Per (h)"],vmin_meta,vmax_meta, periodmap_color, 
       'Metacyc Per (' r'$\mu$' +': ' +str(meta_meanper) +'h)','period (h)', 8, tick = "period")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Period Maps/Metacyc Period Map .tif", format="TIF", dpi=dpi)

#CWT period map
mapper(compiled_period["Avg CWT Per (h)"],vmin_cwt,vmax_cwt, periodmap_color, 'CWT Per ('r'$\mu$' +': ' +str(cwt_meanper) + 'h)',
       'period (h)', 8,tick = "period")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Period Maps/CWT Period Map .tif", format="TIF", dpi=dpi)

#Fourier period map
mapper(compiled_period['FFT Per (h)'],vmin_fft,vmax_fft, periodmap_color, 'FFT Per ('r'$\mu$' +': ' +str(FFT_meanper) +'h)',
 'period (h)', 8, tick = "period")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Period Maps/FFT Period Map .tif", format="TIF", dpi=dpi)

#posslope period map
mapper(compiled_period["XOV + Slope Per (h)"],vmin_xov, vmax_xov, periodmap_color, 
       '+ Slope Per ('r'$\mu$' +': ' +str(pslope_meanper) + 'h)','period (h)', 8, tick = "period")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Period Maps/Pos Slope Period Map .tif", format="TIF", dpi=dpi)

#negslope period map
mapper(compiled_period["XOV - Slope Per (h)"],vmin_xov, vmax_xov, periodmap_color, 
       '- Slope Per ('r'$\mu$' +': ' +str(nslope_meanper) + 'h)','period (h)', 8, tick = "period")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Period Maps/Neg Slope Period Map .tif", format="TIF", dpi=dpi)

#peak period map
mapper(compiled_period["XOV Peak Per (h)"],vmin_xov, vmax_xov, periodmap_color, 
       'Peak Per ('r'$\mu$' +': ' +str(peak_meanper) + 'h)','period (h)', 8, tick = "period")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Period Maps/Peak Period Map .tif", format="TIF", dpi=dpi)

#trough period map
mapper(compiled_period["XOV Trough Per (h)"],vmin_xov, vmax_xov, periodmap_color, 
       'Trough Per ('r'$\mu$' +': ' +str(trough_meanper) + 'h)','period (h)', 8, tick = "period")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Period Maps/Trough Period Map .tif", format="TIF", dpi=dpi)


#metacycle rel phase map
mapper(meta_rel_phase *12/np.pi,-12,12, 
       phasemap_color, 'rPha at ' + str(phasemap_time) + ' h (Metacyc) 'r'$\mu$' +': '  
       + str(meta_meanphaseh) + 'h', 'phase diff (h)', 8, tick = "phase")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/Metacyc rPhase Map .tif", format="TIF", dpi=dpi)

#metacycle abs phase map
mapper(phaset_meta*12/np.pi,0,24, 
       phasemap_color, 'Pha at ' + str(phasemap_time) + 'h (Metacyc) 'r'$\mu$' +': '  
       + str(meta_meanphaseh) + 'h','phase (h)', 8, tick = "phase")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/Metacyc Phase Map .tif", format="TIF", dpi=dpi)

#CWT Rel phase map
mapper(cwt_rel_phase*12/np.pi,-12,12, 
       phasemap_color, 'rPhase at ' +str(phasemap_time) + 'h (CWT) 'r'$\mu$' +': '  
       + str(cwt_meanphaseh) +'h', 'phase diff (h)', 8, tick = "phase")

plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/CWT rPhase Map .tif", format="TIF", dpi=dpi)

#CWT Abs phase map
mapper(phaset_cwt *12/np.pi, 0,24, 
       phasemap_color, 'Phase at ' +str(phasemap_time) + ' h (CWT) ' r'$\mu$' +': '  
       + str(cwt_meanphaseh) +'h', 'phase (h)', 8, tick = "phase")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/CWT Phase Map.tif", format="TIF", dpi=dpi)

#posslope Abs phase map
mapper(posslope_phase *12/np.pi, 0,24, 
       phasemap_color, '+ Slope Phase ' + r'$\mu$' +': '  
       + str(round(posslope_meanphase*12/np.pi,2)) +'h', 'phase (h)', 8, tick = "phase")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/XOV Pos Slope Phase Map.tif", format="TIF", dpi=dpi)

#posslope Rel phase map
mapper(posslope_rel_phase *12/np.pi, 0,24, 
       phasemap_color, '+ Slope rPhase ' + r'$\mu$' +': '  
       + str(round(posslope_meanphase*12/np.pi,2)) +'h', 'phase diff (h)', 8, tick = "phase")

plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/XOV Pos Slope rPhase Map.tif", format="TIF", dpi=dpi)

#negslope Abs phase map
mapper(negslope_phase *12/np.pi, 0,24, 
       phasemap_color, '- Slope Phase ' + r'$\mu$' +': '  
       + str(round(negslope_meanphase*12/np.pi,2)) +'h', 'phase (h)', 8, tick = "phase")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/XOV Neg Slope Phase Map.tif", format="TIF", dpi=dpi)

#negslope Rel phase map
mapper(negslope_rel_phase *12/np.pi, 0,24, 
       phasemap_color, '- Slope rPhase ' + r'$\mu$' +': '  
       + str(round(negslope_meanphase*12/np.pi,2)) +'h', 'phase diff (h)', 8, tick = "phase")

plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/XOV Neg Slope rPhase Map.tif", format="TIF", dpi=dpi)

#Peak Abs phase map
mapper(peak_phase *12/np.pi, 0,24, 
       phasemap_color, 'Peak Phase ' + r'$\mu$' +': '  
       + str(round(peak_meanphase*12/np.pi,2)) +'h', 'phase (h)', 8, tick = "phase")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/XOV Peak Slope Phase Map.tif", format="TIF", dpi=dpi)

#peak Rel phase map
mapper(peak_rel_phase *12/np.pi, 0,24, 
       phasemap_color, 'Peak rPhase ' + r'$\mu$' +': '  
       + str(round(peak_meanphase*12/np.pi,2)) +'h', 'phase (h)', 8, tick = "phase")

plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/XOV Peak rPhase Map.tif", format="TIF", dpi=dpi)


#Troughlope Abs phase map
mapper(trough_phase *12/np.pi, 0,24, 
       phasemap_color, 'Trough Phase ' + r'$\mu$' +': '  
       + str(round(trough_meanphase*12/np.pi,2)) +'h', 'phase (h)', 8, tick = "phase")
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/XOV Trough Phase Map.tif", format="TIF", dpi=dpi)

#posslope Rel phase map
mapper(posslope_rel_phase *12/np.pi, 0,24, 
       phasemap_color, '+ Slope rPhase ' + r'$\mu$' +': '  
       + str(round(trough_meanphase*12/np.pi,2)) +'h', 'phase diff (h)', 8, tick = "phase")

plt.savefig("./Results_" +str(filename_raw)+"/Plots/Maps/Phase Maps/XOV Trough rPhase Map.tif", format="TIF", dpi=dpi)



# Plots period histograms
#reads data
posslope_per = pd.read_csv("./Results_" +str(filename_raw)+"/Crossover pos slope per.csv", index_col=0)
negslope_per = pd.read_csv("./Results_" +str(filename_raw)+"/Crossover neg slope per.csv", index_col=0)
peaks_per = pd.read_csv("./Results_" +str(filename_raw)+"/Crossover peak per.csv", index_col=0)
troughs_per = pd.read_csv("./Results_" +str(filename_raw)+ "/Crossover trough per.csv",  index_col=0)
#average cyclewise period stability
posslope_per_sd = round(posslope_per.mean().std(),2)
negslope_per_sd = round(negslope_per.mean().std(),2)
peaks_per_sd = round(peaks_per.mean().std(),2)
troughs_per_sd = round(troughs_per.mean().std(),2)



if min_per == 'auto' and max_per == 'auto':
    minimum = [min(compiled_period["Metacyc Per (h)"]), min(compiled_period['Avg CWT Per (h)']), 
               min(compiled_period['FFT Per (h)'])]
    minimum = pd.Series(minimum)
    maximum = [max(compiled_period["Metacyc Per (h)"]), max(compiled_period['Avg CWT Per (h)']),
               max(compiled_period['FFT Per (h)'])]
    maximum = pd.Series(maximum)
    axlim_min = min(minimum) - 0.5
    axlim_max = maximum.max()+0.5
else:
    axlim_min = min_per - 0.5
    axlim_max = max_per + 0.5

axticks = np.round(np.linspace(axlim_min, axlim_max, num = 7),1)
axrotation = 45



#correlation plots for period estimation method comparison
pr1 = pearsonr(compiled_period["Metacyc Per (h)"],compiled_period['Avg CWT Per (h)'])
pr1p = round(pr1[1],3)
pr1 = round(pr1[0],2)
pr2 = pearsonr(compiled_period["Metacyc Per (h)"],compiled_period['FFT Per (h)'])
pr2p = round(pr2[1],3)
pr2 = round(pr2[0],2)
pr3 = pearsonr(compiled_period['Avg CWT Per (h)'],compiled_period['FFT Per (h)'])
pr3p = round(pr3[1],3)
pr3 = round(pr3[0],2)

fig, axs = plt.subplots(1, 3, sharey=False, tight_layout=True, figsize=(8,2.75))
plt.setp(axs, xticks =(axticks), yticks = (axticks), xlim = (axlim_min, axlim_max), ylim = (axlim_min, axlim_max))

axs[0].scatter(compiled_period["Metacyc Per (h)"],compiled_period['FFT Per (h)'],color = "tab:green", s=2)
axs[0].set_title("corr = " + str(pr2) + "; p = " + str(pr2p))
axs[0].set_ylabel('FFT period (h)')
axs[0].set_xlabel('Metacyc period (h)') 
axs[0].set_xticklabels(axticks, rotation = axrotation)
axs[1].scatter(compiled_period["Metacyc Per (h)"],compiled_period['Avg CWT Per (h)'], color = "tab:red", s = 2)
axs[1].set_title("corr = " + str(pr1) + "; p = " + str(pr1p))
axs[1].set_ylabel('CWT period (h)')
axs[1].set_xlabel('Metacyc period (h)')
axs[1].set_xticklabels(axticks, rotation = axrotation)
axs[2].scatter(compiled_period['Avg CWT Per (h)'],compiled_period['FFT Per (h)'], color = "tab:blue", s=2) 
axs[2].set_title("corr = " + str(pr3) + "; p = " + str(pr3p))
axs[2].set_xticklabels(axticks, rotation = axrotation)
axs[2].set_xlabel('CWT period(h)')
axs[2].set_ylabel('FFT period (h)')
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Correlation Plots 1.tif", format="TIF", dpi=dpi)

temp1 = pd.concat([compiled_period["Avg CWT Per (h)"],compiled_period['XOV + Slope Per (h)']], axis = 1).dropna()
pr4 = pearsonr(temp1["Avg CWT Per (h)"],temp1['XOV + Slope Per (h)'])
pr4p = round(pr4[1],3)
pr4 = round(pr4[0],2)
temp2 = pd.concat([compiled_period["Avg CWT Per (h)"],compiled_period['XOV - Slope Per (h)']], axis = 1).dropna()
pr5 = pearsonr(temp2["Avg CWT Per (h)"],temp2['XOV - Slope Per (h)'])
pr5p = round(pr5[1],3)
pr5 = round(pr5[0],2)
temp3 = pd.concat([compiled_period["Avg CWT Per (h)"],compiled_period['XOV Peak Per (h)']], axis = 1).dropna()
pr6 = pearsonr(temp3["Avg CWT Per (h)"],temp3['XOV Peak Per (h)'])
pr6p = round(pr6[1],3)
pr6 = round(pr6[0],2)
temp4 = pd.concat([compiled_period["Avg CWT Per (h)"],
                          compiled_period['XOV Trough Per (h)']], axis = 1).dropna()
pr7 = pearsonr(temp4["Avg CWT Per (h)"],temp4['XOV Trough Per (h)'])
pr7p = round(pr7[1],3)
pr7 = round(pr7[0],2)

fig, axs = plt.subplots(1, 4, sharey=False, tight_layout=True, figsize=(12,3))
plt.setp(axs, xticks =(axticks), yticks = (axticks), xlim = (axlim_min, axlim_max), ylim = (axlim_min, axlim_max))
axs[0].scatter(compiled_period["Avg CWT Per (h)"],compiled_period['XOV + Slope Per (h)'],  color = "tab:orange", s=2) 
axs[0].set_title("corr = " + str(pr4) + "; p = " + str(pr4p))
axs[0].set_xticklabels(axticks, rotation = axrotation)
axs[0].set_xlabel('CWT period (h)')
axs[0].set_ylabel('XOV + Slope CWT Per (h)')
axs[1].scatter(compiled_period["Avg CWT Per (h)"],compiled_period['XOV - Slope Per (h)'], color = "tab:purple", s=2) 
axs[1].set_title("corr = " + str(pr5) + "; p = " + str(pr5p))
axs[1].set_xticklabels(axticks, rotation = axrotation)
axs[1].set_xlabel('CWT period (h)')
axs[1].set_ylabel('XOV - Slope Per (h)')
axs[2].scatter(compiled_period["Avg CWT Per (h)"],compiled_period['XOV Peak Per (h)'], color = "tab:brown", s=2) 
axs[2].set_title("corr = " + str(pr6) + "; p = " + str(pr6p))
axs[2].set_xticklabels(axticks, rotation = axrotation)
axs[2].set_xlabel('CWT period (h)')
axs[2].set_ylabel('XOV Peak Per (h)')
axs[3].scatter(compiled_period["Avg CWT Per (h)"],
                          compiled_period['XOV Trough Per (h)'], color = "tab:brown", s=2) 
axs[3].set_title("corr = " + str(pr7) + "; p = " + str(pr7p))
axs[3].set_xticklabels(axticks, rotation = axrotation)
axs[3].set_xlabel('CWT period (h)')
axs[3].set_ylabel('XOV Trough Per (h)')

plt.savefig("./Results_" +str(filename_raw)+"/Plots/Correlation Plots 2.tif", format="TIF", dpi=dpi)


#histograms of periods
fig, axs = plt.subplots(2,4 , sharey=False, tight_layout=True, figsize=(12,6))
axs[0,0].hist(compiled_period["Metacyc Per (h)"],color = "tab:red", alpha = 0.9)
axs[0,0].set_title(label = 'Metacyc Per ' r'$\mu$' +': ' + 
                   str(meta_meanper) + 'h', fontsize = '10')
axs[0,0].set_ylabel('# of cells')
axs[0,0].set_xlabel('period (h)')
axs[0,0].set_xlim(axlim_min, axlim_max)
axs[0,0].set_xticks(axticks)
axs[0,0].set_xticklabels(axticks, rotation = axrotation)

axs[0,1].hist(compiled_period["FFT Per (h)"], color = "tab:red", alpha = 0.8) 
axs[0,1].set_title(label = 'FFT Per 'r'$\mu$' +':' + 
                   str(FFT_meanper) + 'h', fontsize = '10')
axs[0,1].set_ylabel('# of cells')
axs[0,1].set_xlim(axlim_min, axlim_max)
axs[0,1].set_xticks(axticks)
axs[0,1].set_xticklabels(axticks, rotation = axrotation)
axs[0,1].set_xlabel('period (h)')

axs[0,2].hist(compiled_period["Avg CWT Per (h)"], color = "tab:red", alpha = 0.9)
axs[0,2].set_ylabel('# of cells')
axs[0,2].set_xlim(axlim_min, axlim_max)
axs[0,2].set_xticks(axticks)
axs[0,2].set_xticklabels(axticks, rotation = axrotation)
axs[0,2].set_title(label = 'CWT Per ' r'$\mu$' +': ' +str(cwt_meanper) + 'h', fontsize = '10')
axs[0,2].set_xlabel('period (h)')

axs[0,3].hist(compiled_period['XOV + Slope Per (h)'], color = "tab:red", alpha = 0.9) 
axs[0,3].set_title(label = '+ slope Per 'r'$\mu$' +': ' +
                   str(pslope_meanper) + 'h', fontsize = '10')
axs[0,3].set_ylabel('# of cells')
axs[0,3].set_xlim(axlim_min, axlim_max)
axs[0,3].set_xticks(axticks)
axs[0,3].set_xticklabels(axticks, rotation = axrotation)
axs[0,3].set_xlabel('period (h)')

axs[1,0].hist(compiled_period['XOV - Slope Per (h)'],  color = "tab:red", alpha = 0.9) 
axs[1,0].set_title(label = '- slope Per 'r'$\mu$' +': ' +
                   str(nslope_meanper) + 'h', fontsize = '10')
axs[1,0].set_ylabel('# of cells')
axs[1,0].set_xlim(axlim_min, axlim_max)
axs[1,0].set_xticks(axticks)
axs[1,0].set_xticklabels(axticks, rotation = axrotation)
axs[1,0].set_xlabel('period (h)')

axs[1,1].hist(compiled_period['XOV Peak Per (h)'],  color = "tab:red", alpha = 0.9) 
axs[1,1].set_title(label = 'Peak-Peak Per 'r'$\mu$' +': ' +
                   str(peak_meanper) + 'h', fontsize = '10')
axs[1,1].set_ylabel('# of cells')
axs[1,1].set_xlim(axlim_min, axlim_max)
axs[1,1].set_xticks(axticks)
axs[1,1].set_xticklabels(axticks, rotation = axrotation)
axs[1,1].set_xlabel('period (h)')


axs[1,2].hist(compiled_period['XOV Trough Per (h)'],  color = "tab:red", alpha = 0.9) 
axs[1,2].set_title(label = 'Trough-Trough Per 'r'$\mu$' +': ' +
                   str(trough_meanper) + 'h', fontsize = '10')
axs[1,2].set_ylabel('# of cells')
axs[1,2].set_xlim(axlim_min, axlim_max)
axs[1,2].set_xticks(axticks)
axs[1,2].set_xticklabels(axticks, rotation = axrotation)
axs[1,2].set_xlabel('period (h)')

axs[1,3].bar([1,2,3,4], [posslope_per_sd, negslope_per_sd, peaks_per_sd, troughs_per_sd], color = "tab:blue", alpha = 0.9) 
axs[1,3].set_title(label = 'Inter-cycle Per SD', fontsize = '10')
axs[1,3].set_ylabel('Standard Deviation')
axs[1,3].set_xticks([1,2,3,4])
axs[1,3].set_xticklabels(['+ slope', '- slope', 'Peak', 'Trough'], rotation = axrotation)
axs[1,3].set_xlabel('Phase Marker')
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Period Histograms.tif", format="TIF", dpi=dpi)

# Plots Amp histograms
fig, axs = plt.subplots(1,3 , sharey=False, tight_layout=True, figsize=(7,3))
axs[0].hist(compiled_amplitude['Metacyc Amp'], color = "tab:green", alpha = 0.8) 
axs[0].set_title(label = 'Metacyc Amp '+ r'$\mu$' +': ' +str(meta_meanamp), fontsize = '8')
axs[0].set_xlabel('Amplitude')
axs[0].set_ylabel('# of cells')
axs[0].tick_params(axis='x', labelrotation=axrotation)
axs[1].hist(compiled_amplitude['Avg CWT Amp'],  color = "tab:green", alpha = 0.8) 
axs[1].set_title(label = 'Avg CWT Amp '+ r'$\mu$' 
                   +': ' +str(cwt_meanamp), fontsize = '8')
axs[1].set_xlabel('Amplitude')
axs[1].set_ylabel('# of cells')
axs[1].tick_params(axis='x', labelrotation=axrotation)
axs[2].hist(compiled_amplitude['XOV Amp'], color = "tab:green", alpha = 0.8) 
axs[2].set_title(label = 'Crossover Amp '+ r'$\mu$' 
                   +': ' +str(xov_meanamp), fontsize = '8')
axs[2].set_xlabel('Amplitude')
axs[2].set_ylabel('# of cells')
axs[2].tick_params(axis='x', labelrotation=axrotation)
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Amplitude Histograms.tif", format="TIF", dpi=dpi)

# Plots Phase histograms
fig, axs = plt.subplots(2,6 , sharey=False, tight_layout=True, figsize=(12,4))
axs[0,0].hist(round(phaset_meta*12/np.pi,2),  color = "tab:blue", alpha = 0.9)
axs[0,0].set_title(label = 'Metacyc Pha at ' +str(phasemap_time) + "h "r'$\mu$' +': ' +str(round(meta_meanphaseh,2)) + 'h', 
                   fontsize = '8')
axs[0,0].set_xlabel('Phase (h)')
axs[0,0].set_ylabel('# of cells')
axs[0,0].set_xlim(0, 24)
axs[0,0].set_xticks([0,4,8,12,16,20,24])
axs[0,0].set_xticklabels([0,4,8,12,16,20,24], rotation = axrotation)
axs[1,0].hist(round(meta_rel_phase*12/np.pi,2), color = "tab:blue", alpha = 0.9)
axs[1,0].set_title(label = 'Metacyc rPha at ' +str(phasemap_time) + "h "r'$\mu$' +': ' +str(round(meta_meanphaseh,2)) + 'h', 
                   fontsize = '8')
axs[1,0].set_xlabel('Phase Diff (h)')
axs[1,0].set_ylabel('# of cells')
axs[1,0].set_xlim(-12, 12)
axs[1,0].set_xticks([-12,-8,-4,0,4,8,12])
axs[1,0].set_xticklabels([-12,-8,-4,0,4,8,12], rotation = axrotation)
axs[0,1].hist(round(phaset_cwt*12/np.pi,2), color = "tab:blue", alpha = 0.9) 
axs[0,1].set_xlim(0, 24)
axs[0,1].set_xticks([0,4,8,12,16,20,24])
axs[0,1].set_xticklabels([0,4,8,12,16,20,24], rotation = axrotation)
axs[0,1].set_title(label = 'CWT Pha at ' +str(phasemap_time) + "h "r'$\mu$' +': ' +
                   str(round(cwt_meanphaseh,2)) + 'h', fontsize = '8')
axs[0,1].set_xlabel('Phase (h)')
axs[0,1].set_ylabel('# of cells')
axs[1,1].hist(round(cwt_rel_phase*12/np.pi,2), color = "tab:blue", alpha = 0.9) 
axs[1,1].set_xlim(-12, 12)
axs[1,1].set_xticks([-12,-8,-4,0,4,8,12])
axs[1,1].set_title(label = 'CWT rPha at ' +str(phasemap_time) + "h "r'$\mu$' +': ' +
                   str(round(cwt_meanphaseh,2)) + 'h', fontsize = '8')
axs[1,1].set_xlabel('Phase (h)')
axs[1,1].set_ylabel('# of cells')
axs[1,1].set_xticklabels([-12,-8,-4,0,4,8,12], rotation = axrotation)

axs[0,2].hist(xov_phase['XOV + Slope Pha (h)'], color = "tab:blue", alpha = 0.9) 
axs[0,2].set_xlim(0, 24)
axs[0,2].set_xticks([0,4,8,12,16,20,24])
axs[0,2].set_xticklabels([0,4,8,12,16,20,24], rotation = axrotation)
axs[0,2].set_title(label = '+ Slope Pha ' + r'$\mu$' +': ' +
                   str(round(posslope_meanphase*12/np.pi,2)) + 'h', fontsize = '8')
axs[0,2].set_xlabel('Phase (h)')
axs[0,2].set_ylabel('# of cells')

axs[1,2].hist(xov_phase['XOV + Slope rPha (h)'], color = "tab:blue", alpha = 0.9) 
axs[1,2].set_xlim(-12, 12)
axs[1,2].set_xticks([-12,-8,-4,0,4,8,12])
axs[1,2].set_title(label = '+ Slope rPha ' + r'$\mu$' +': ' +
                   str(round(posslope_meanphase*12/np.pi,2)) + 'h', fontsize = '8')
axs[1,2].set_xlabel('Phase Diff (h)')
axs[1,2].set_ylabel('# of cells')
axs[1,2].set_xticklabels([-12,-8,-4,0,4,8,12], rotation = axrotation)

axs[0,3].hist(xov_phase['XOV - Slope Pha (h)'], color = "tab:blue", alpha = 0.9) 
axs[0,3].set_xlim(0, 24)
axs[0,3].set_xticks([0,4,8,12,16,20,24])
axs[0,3].set_xticklabels([0,4,8,12,16,20,24], rotation = axrotation)
axs[0,3].set_title(label = '- Slope Pha ' + r'$\mu$' +': ' +
                   str(round(negslope_meanphase*12/np.pi,2)) + 'h', fontsize = '8')
axs[0,3].set_xlabel('Phase (h)')
axs[0,3].set_ylabel('# of cells')

axs[1,3].hist(xov_phase['XOV - Slope rPha (h)'], color = "tab:blue", alpha = 0.9) 
axs[1,3].set_xlim(-12, 12)
axs[1,3].set_xticks([-12,-8,-4,0,4,8,12])
axs[1,3].set_xticklabels([-12,-8,-4,0,4,8,12], rotation = axrotation)
axs[1,3].set_title(label = '- Slope rPha ' + r'$\mu$' +': ' +
                   str(round(negslope_meanphase*12/np.pi,2)) + 'h', fontsize = '8')
axs[1,3].set_xlabel('Phase Diff (h)')
axs[1,3].set_ylabel('# of cells')


axs[0,4].hist(xov_phase['XOV Peak Pha (h)'], color = "tab:blue", alpha = 0.9) 
axs[0,4].set_xlim(0, 24)
axs[0,4].set_xticks([0,4,8,12,16,20,24])
axs[0,4].set_xticklabels([0,4,8,12,16,20,24], rotation = axrotation)
axs[0,4].set_title(label = 'Peak Pha ' + r'$\mu$' +': ' +
                   str(round(peak_meanphase*12/np.pi,2)) + 'h', fontsize = '8')
axs[0,4].set_xlabel('Phase (h)')
axs[0,4].set_ylabel('# of cells')

axs[1,4].hist(xov_phase['XOV Peak rPha (h)'], color = "tab:blue", alpha = 0.9) 
axs[1,4].set_xlim(-12, 12)
axs[1,4].set_xticks([-12,-8,-4,0,4,8,12])
axs[1,4].set_xticklabels([-12,-8,-4,0,4,8,12], rotation = axrotation)
axs[1,4].set_title(label = 'Peak rPha ' + r'$\mu$' +': ' +
                   str(round(peak_meanphase*12/np.pi,2)) + 'h', fontsize = '8')
axs[1,4].set_xlabel('Phase Diff (h)')
axs[1,4].set_ylabel('# of cells')

axs[0,5].hist(xov_phase['XOV Trough Pha (h)'], color = "tab:blue", alpha = 0.9) 
axs[0,5].set_xlim(0, 24)
axs[0,5].set_xticks([0,4,8,12,16,20,24])
axs[0,5].set_xticklabels([0,4,8,12,16,20,24], rotation = axrotation)
axs[0,5].set_title(label = 'Trough Pha ' + r'$\mu$' +': ' +
                   str(round(trough_meanphase*12/np.pi,2)) + 'h', fontsize = '8')
axs[0,5].set_xlabel('Phase (h)')
axs[0,5].set_ylabel('# of cells')

axs[1,5].hist(xov_phase['XOV Trough rPha (h)'], color = "tab:blue", alpha = 0.9) 
axs[1,5].set_xlim(-12, 12)
axs[1,5].set_xticks([-12,-8,-4,0,4,8,12])
axs[1,5].set_xticklabels([-12,-8,-4,0,4,8,12], rotation = axrotation)
axs[1,5].set_title(label = 'Trough rPha ' + r'$\mu$' +': ' +
                   str(round(trough_meanphase*12/np.pi,2)) + 'h', fontsize = '8')
axs[1,5].set_xlabel('Phase Diff (h)')
axs[1,5].set_ylabel('# of cells')
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Phase Histograms.tif", format="TIF", dpi=dpi)




#Plots cyclewise crossover and cwt amp  data
#calculates mean across cells for every cycle
xov_cycwise_meanamp = np.array(xov_amp_cycwise.mean())
xov_cycwise_meanamp=xov_cycwise_meanamp[~np.isnan(xov_cycwise_meanamp)]
xov_cycwise_ampstd = np.array(xov_amp_cycwise.std()) 
xov_cycwise_ampstd= xov_cycwise_ampstd[~np.isnan(xov_cycwise_ampstd)]
xov_cycwise_normamp = np.round(100*(xov_cycwise_meanamp/xov_cycwise_meanamp[0]),2)
cwt_amplitude = pd.DataFrame(amp_envelope).T
ncycles_cwt = int(len(cwt_amplitude)*sampling_interval_mins/(24*60))
cwt_time = np.arange(0, int((24*ncycles_com*(sampling_interval_mins/60))+1), int(24*sampling_interval_mins/60))
legend =[]
for i in np.arange(1, ncycles_com+1):
    legend.append('day ' + str(i))

daylength = int(24*60/sampling_interval_mins)
daylength = daylength*np.arange(0,ncycles_cwt+1)
cwt_cycwise_amp=[]
for i in range(len(daylength)-1):
    cwt_cycwise_amp.append(cwt_amplitude.loc[daylength[i]:daylength[i+1]].mean().mean())
cwt_cycwise_amp = pd.DataFrame(cwt_cycwise_amp)
cwt_cycwise_norm_amp = np.round(100*(cwt_cycwise_amp[0]/cwt_cycwise_amp[0][0]),2)


fig, axs = plt.subplots(2,2 , sharey=False, tight_layout=True, figsize=(6,5))
axs[0,0].bar(np.arange(1, len(cwt_cycwise_amp)+1),cwt_cycwise_amp[0] ,color = 'tab:orange')
#axs[0,1].set_ylim(0,max(xov_ampstd)+max(xov_meanamp)+)
axs[0,0].set_xlabel("Day", fontsize = '10')
axs[0,0].set_ylabel("Amp", fontsize = '10')
axs[0,0].set_title('Daywise Amp (CWT)', fontsize = '10')
axs[0,1].bar(np.arange(1, len(pd.DataFrame(cwt_cycwise_norm_amp)[0])+1),pd.DataFrame(cwt_cycwise_norm_amp)[0],color = 'tab:orange')
#axs[0,1].set_ylim(0,max(xov_ampstd)+max(xov_meanamp)+)
axs[0,1].set_xlabel("Day", fontsize = '10')
axs[0,1].set_ylabel("Normalized Amp", fontsize = '10')
axs[0,1].set_title('Daywise Norm Amp (CWT)', fontsize = '10')
axs[1,0].bar(np.arange(1, len(xov_cycwise_meanamp)+1),xov_cycwise_meanamp,color = 'tab:blue')
#axs[0,1].set_ylim(0,max(xov_ampstd)+max(xov_meanamp)+)
axs[1,0].set_xlabel("Cycle", fontsize = '10')
axs[1,0].set_ylabel("Peak-Trough Amp", fontsize = '10')
axs[1,0].set_title('Cycwise Amp (crossover)', fontsize = '10')
axs[1,1].bar(np.arange(1, len(xov_cycwise_meanamp)+1), xov_cycwise_normamp  ,color = 'tab:blue')
axs[1,1].set_xlabel("Cycle", fontsize = '10')
axs[1,1].set_ylabel("Normalized Peak-Trough Amp", fontsize = '10')
axs[1,1].set_title('Cycwise Norm Amp(crossover)', fontsize = '10')

plt.savefig("./Results_" +str(filename_raw)+"/Plots/Daywise Amplitude.tif", format="TIF", dpi =dpi)
plt.show()


#Plots Ensemble dynamics
fig,axs = plt.subplots(2,2, figsize = (8, 5), tight_layout = True)
xaxs = round(pd.DataFrame(time)/60,1).T
xaxs = xaxs.squeeze(axis = 1)
#plt.setp(axs, xlim = (min(xaxs)-1, max(xaxs)+1))
axs[0,0].plot(xaxs,avg_ts, color = 'tab:blue')
axs[0,0].plot(xaxs,avg_inst_amp, color = 'tab:orange')
axs[0,0].set_xlabel("Time(h)", fontsize = '10')
axs[0,0].set_ylabel("signal intensity", fontsize = '10')
axs[0,0].set_title('Avg Signal and Amplitude', fontsize = '13')
axs[0,0].legend(['Avg Time Series', 'Avg Inst amplitude'], prop={'size': 4})
axs[0,0].vlines([(60/sampling_interval_mins)*6, xaxs.loc[len(xaxs)-((60/sampling_interval_mins)*6)]],
                min(avg_ts) , max(avg_ts), colors = 'grey', linestyles = 'dashed')
axs[0,1].plot(xaxs,avg_inst_per, color = 'tab:green')
axs[0,1].vlines([(60/sampling_interval_mins)*6, xaxs.loc[len(xaxs)-((60/sampling_interval_mins)*6)]],
                0,30, colors = 'grey', linestyles = 'dashed')
axs[0,1].set_xlabel("Time(h)", fontsize = '10')
axs[0,1].set_ylabel("Period (h)", fontsize = '10')
axs[0,1].set_ylim(min_per, max_per)
axs[0,1].set_title('Avg Instantaneous Period', fontsize = '13')
axs[1,0].plot(xaxs, r1, color = 'tab:brown')
axs[1,0].set_ylim(0,1.2)
axs[1,0].set_xlabel("Time(h)", fontsize = '10')
axs[1,0].set_ylabel("Order Parameter", fontsize = '10')
axs[1,0].set_title('First Order SI', fontsize = '13')
axs[1,0].vlines([(60/sampling_interval_mins)*6, xaxs.loc[len(xaxs)-((60/sampling_interval_mins)*6)]],
                0,1.2, colors = 'grey', linestyles = 'dashed')
axs[1,1].plot(xaxs, r2, color = 'tab:brown')
axs[1,1].set_ylim(0,1.2)
axs[1,1].set_xlabel("Time(h)", fontsize = '10')
axs[1,1].set_ylabel("Order Parameter",fontsize = '10')
axs[1,1].set_title('Second Order SI', fontsize = '13')
axs[1,1].vlines([(60/sampling_interval_mins)*6, xaxs.loc[len(xaxs)-((60/sampling_interval_mins)*6)]],
                0,1.2, colors = 'grey', linestyles = 'dashed')
plt.savefig("./Results_" +str(filename_raw)+"/Plots/Ensemble Dynamics.tif", format="TIF", dpi =dpi)
plt.show()

#Phase evolution plots

#maps
for i in range(len(time_range)):
    z = cwt_phase[i] * 12/np.pi
    mapper(z,0,24, phasemap_color, "CWT Phase (t =" + str(round(time_range[0][i]/60,2)) + "h)",
           "phase (h)", 6, tick = "phase")
    plt.savefig("./Results_" +str(filename_raw)+"/Plots/Stack Plots/CWT Abs Phase Evol Map/Phase(t=" 
                + str(round(time_range[0][i]/60,2)) +").tif", format="TIF", dpi=dpi)
    plt.close('all')
    #print("saving phase(time" + str(i/60))

#polar plots (Abs Phase)
for i in range(len(time_range)):
    theta1 = cwt_phase[i]
    #theta_mean = np.arctan2(np.sin(theta1).sum(),np.cos(theta1).sum())
    #theta = (theta1 - theta_mean)
    r = np.repeat(1, len(theta1))
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(theta1,r,'o', color = 'tab:green') 
    ax.set_theta_zero_location("N")
    ax.axes.yaxis.set_visible(False)
    ax.grid(True)
    ax.set_theta_direction(-1)
    ax.set_title("CWT Phase (t = " +str(round(time_range[0][i]/60,2)) +" h)", va='bottom')
    plt.savefig("./Results_" +str(filename_raw)+"/Plots/Stack Plots/CWT Abs Phase Evol Polar/Phase(t="
                + str(round(time_range[0][i]/60,2)) +").tif", format="TIF", dpi=dpi)
    plt.close('all')

#polar plots (Rel Phase)
for i in range(len(time_range)):
    theta1 = cwt_phase[i]
    theta_mean = np.arctan2(np.sin(theta1).sum(),np.cos(theta1).sum())
    theta = (theta1 - theta_mean)
    r = np.repeat(1, len(theta1))
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(theta,r,'o', color = 'tab:green') 
    ax.set_theta_zero_location("N")
    ax.axes.yaxis.set_visible(False)
    ax.grid(True)
    ax.set_theta_direction(-1)
    ax.set_title("CWT Rel Phase (t = " +str(round(time_range[0][i]/60,2)) +" h)", va='bottom')
    plt.savefig("./Results_" +str(filename_raw)+"/Plots/Stack Plots/CWT Rel Phase Evol Polar/Phase(t="
                + str(round(time_range[0][i]/60,2)) +").tif", format="TIF", dpi=dpi)
    plt.close('all')

#polar histograms (CWT abs Phase)
for i in range(len(time_range)):
    degrees = np.rad2deg(cwt_phase[i])
    radians = cwt_phase[i]
    bin_size = 5
    a , b=np.histogram(degrees, bins=np.arange(0, 360+bin_size, bin_size))
    centers = np.deg2rad(np.ediff1d(b)//2 + b[:-1])
    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111, projection='polar')
    ax.bar(centers, a, width=np.deg2rad(bin_size), bottom=0.0, color='tab:orange',alpha = 0.9)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.axes.yaxis.set_visible(False)
    ax.axes.xaxis.set_visible(True)
    ax.set_title("CWT Phase (t = " +str(round(time_range[0][i]/60,2)) +" h)", va='bottom')
    plt.savefig("./Results_" +str(filename_raw)+"/Plots/Stack Plots/Polar Histograms/CWT Abs Phase Polar Histograms/Phase Dist(t="+ str(round(time_range[0][i]/60,2)) +").tif", format="TIF", dpi=dpi)
    plt.close('all')
 
    
#Plots heatmap of CWT phases. 'method' and 'metric' can be changed based on user requirements for better plots
#https://docs.scipy.org/doc/scipy/reference/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
#https://docs.scipy.org/doc/scipy/reference/reference/generated/scipy.spatial.distance.pdist.html#scipy.spatial.distance.pdist
sns.clustermap(pd.DataFrame(xov_smoothed), method='average', metric='correlation',
               z_score = 0, row_cluster=True, col_cluster=False, row_linkage=None, 
               cmap = clustermap_color, cbar_pos = None)
plt.savefig("./Results_" +str(filename_raw)+ "/Plots/HeatMap.tif", format = "TIF", dpi=dpi)

if os.path.exists("./Results_" +str(filename_raw) + "/metadata.csv"):
    os.remove("./Results_" +str(filename_raw) + "/metadata.csv")
else:
    shutil.move('./metadata.csv', "./Results_" +str(filename_raw))
if os.path.exists("./Results_" +str(filename_raw) +"/metaout"):
    shutil.rmtree("./Results_" +str(filename_raw) +"/metaout")
shutil.move("./metaout", "./Results_" +str(filename_raw)) 
cache = [sampling_interval_mins, min_per, max_per, phasemap_time]
cache = pd.DataFrame(cache).T
cache.columns = ['sampling_interval_mins', 'min_per', 'max_per', 'phasemap_time']
cache.to_csv("./Results_" +str(filename_raw)+"/cache/cache.csv")

