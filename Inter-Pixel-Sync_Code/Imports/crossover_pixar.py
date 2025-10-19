# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 18:09:15 2021

@author: Nikhil
"""
#Crossover analysis
#1. takes detrended timeseries as xov_df
#2. calculates 24 h rolling woindow average as xov_roll_avg
#3. takes difference of xov_df and xov_roll_avg
#4. looks for change in sign which indicates crossover
#5. Calculates avg time period using both crossings (upward to downward crossover and other way around)
from scipy.signal import find_peaks
image = io.imread(filename_full)
dim = rescale(image, image_scaling_factor, anti_aliasing=True)
im = np.empty((image.shape[0],dim.shape[1],dim.shape[2]))
for i in range(image.shape[0]):
    im [i] = rescale(image[i], image_scaling_factor, anti_aliasing=True)
im = gaussian(im, sigma=3) #- gaussian(im, sigma=20.0)

#reads gaussian filtered (smoothed) values for pixels reported as rhythmic by metacycle
pix_val_ts = []
for i in range(len(rhythmic_pixs)):
    pix_val_ts.append(im[:,rhythmic_pixs['Y'].loc[i],rhythmic_pixs['X'].loc[i]])
pix_val_ts = pd.DataFrame(pix_val_ts).T/pd.DataFrame(pix_val_ts).mean(axis = 1)

periods = np.round(np.linspace(min_per*60, (max_per)*60, nper), 2)
xov_smoothed = []

wAn = WAnalyzer(periods, sampling_interval_mins, time_unit_label = "min")
for i in range(len(pix_val_ts.columns)):
    filtered = wAn.sinc_detrend(signal = pix_val_ts[i], T_c = T_c*60)
    xov_smoothed.append(filtered)


xov_df = np.array(xov_smoothed)
pad = int((60/sampling_interval_mins)*12)
roll_win = int(pad*2)
xov_df= pd.DataFrame(np.pad(xov_df, ((0,0),(pad,pad)), 'constant')).reset_index(drop = True).T
xov_roll_avg = xov_df.rolling(roll_win, center = True).mean()
xov_roll_avg = xov_roll_avg[int(roll_win/2):
                            int(len(xov_roll_avg) - int(roll_win/2))].reset_index(drop = True)
xov_df = xov_df[int(roll_win/2):int(len(xov_df) - int(roll_win/2))].reset_index(drop = True)

diff= xov_df - xov_roll_avg
sign = np.sign(diff)
diff.insert(0, 'Time',np.arange(0, len(xov_df)*sampling_interval_mins, sampling_interval_mins), True)
posslope_period = []
negslope_period = []
negslope_phase = []
posslope_phase = []
for k in range(len(sign.columns)):
    negslope_xov = []
    posslope_xov = []
    negslope_per = []
    posslope_per = []

    for i in range(len(sign[k])-1):
        if sign.iloc[i,k] == 1.0 and sign.iloc[i+1,k] == -1.0:
            negslope_xov.append(diff['Time'].loc[i])
        elif sign.iloc[i,k] == -1.0 and sign.iloc[i+1,k] == 1.0:
            posslope_xov.append(diff['Time'].loc[i])
        else:
            continue
    negslope_phase.append(negslope_xov)
    posslope_phase.append(posslope_xov)
    for j in range(len(negslope_xov)-1):
            negslope_per.append(negslope_xov[j+1] - negslope_xov[j])
    negslope_period.append(negslope_per)
            
    for j in range(len(posslope_xov)-1):
            posslope_per.append(posslope_xov[j+1] - posslope_xov[j])
    posslope_period.append(posslope_per)

cycwise_negslope_phase = pd.DataFrame(negslope_phase) /60
cycwise_posslope_phase = pd.DataFrame(posslope_phase) /60

avg_negslope_phase = round(cycwise_negslope_phase.mean(axis = 1) %24, 2)
avg_posslope_phase = round(cycwise_posslope_phase.mean(axis = 1) %24,2)

#calculates peaks and trough periods and peak-trough ampplitude
peaks_period =[]
troughs_period = []
peak_phase = []
trough_phase = []
peaks_val = pd.DataFrame()
troughs_val = pd.DataFrame()
for i in range(len(xov_df.columns)):
    peak_per_idx = []
    trough_per_idx = []
    peak_per = []
    trough_per = []
    peak_sig = []
    trough_sig = []
    signal = xov_df[i]
    peak_idx, _ = find_peaks(signal)
    trough_idx, _ = find_peaks(-signal)
    peak_phase.append(peak_idx)
    trough_phase.append(trough_idx)
    peak_per_idx.append(diff['Time'].loc[peak_idx].reset_index(drop = True))
    trough_per_idx.append(diff['Time'].loc[trough_idx].reset_index(drop = True))
    peak_sig.append(signal[peak_idx].reset_index(drop = True))
    trough_sig.append(signal[trough_idx].reset_index(drop = True))
    peaks_val = pd.concat([peaks_val, pd.DataFrame(peak_sig)])
    troughs_val = pd.concat([troughs_val, pd.DataFrame(trough_sig)])
    for j in range(len(peak_per_idx[0])-1):
        peak_per.append(peak_per_idx[0][j+1] - peak_per_idx[0][j]) 
    peaks_period.append(peak_per)
    
    for j in range(len(trough_per_idx[0])-1):
        trough_per.append(trough_per_idx[0][j+1] - trough_per_idx[0][j])
    troughs_period.append(trough_per)

cycwise_peak_phase = pd.DataFrame(peak_phase)*sampling_interval_mins/60
cycwise_trough_phase = pd.DataFrame(trough_phase)*sampling_interval_mins/60
avg_peak_phase = round(cycwise_peak_phase.mean(axis = 1)%24,2)
avg_trough_phase = round(cycwise_trough_phase.mean(axis = 1)%24,2) 

    
#converts periods from minutes to hours and concatenates periods estimated by four methods
xov_posslope_period = round(pd.DataFrame(posslope_period)/60., 4)
xov_posslope_period[xov_posslope_period> max_per] = np.nan
xov_posslope_period[xov_posslope_period< min_per] = np.nan

xov_negslope_period = round(pd.DataFrame(negslope_period)/60.,4)
xov_negslope_period[xov_negslope_period> max_per] = np.nan
xov_negslope_period[xov_negslope_period < min_per] = np.nan
  
peaks_period = round(pd.DataFrame(peaks_period)/60., 4)
peaks_period [peaks_period> max_per] = np.nan
peaks_period [peaks_period< min_per] = np.nan
troughs_period = round(pd.DataFrame(troughs_period)/60., 4)    
troughs_period[troughs_period > max_per] = np.nan
troughs_period[troughs_period < min_per] = np.nan




#calculates relative phases at time t

peak_phase = avg_peak_phase*np.pi/12               
peak_meanphase = (np.arctan2(np.sin(peak_phase).sum(),np.cos(peak_phase).sum()))%(2*np.pi)
peak_rel_phase = peak_phase - peak_meanphase

for i in range(len(peak_rel_phase)):
    if peak_rel_phase[i] > np.pi:
        peak_rel_phase[i] = peak_rel_phase[i]%-np.pi
    elif peak_rel_phase[i] < -1* np.pi:
        peak_rel_phase[i] =  peak_rel_phase[i]%np.pi
    else: continue

trough_phase = avg_trough_phase*np.pi/12               
trough_meanphase = (np.arctan2(np.sin(trough_phase).sum(),np.cos(trough_phase).sum()))%(2*np.pi)
trough_rel_phase = trough_phase - trough_meanphase

for i in range(len(trough_rel_phase)):
    if trough_rel_phase[i] > np.pi:
        trough_rel_phase[i] = trough_rel_phase[i]%-np.pi
    elif trough_rel_phase[i] < -1* np.pi:
        trough_rel_phase[i] =  trough_rel_phase[i]%np.pi
    else: continue

posslope_phase = avg_posslope_phase*np.pi/12               
posslope_meanphase = (np.arctan2(np.sin(posslope_phase).sum(),np.cos(posslope_phase).sum()))%(2*np.pi)
posslope_rel_phase = posslope_phase - posslope_meanphase

for i in range(len(posslope_rel_phase)):
    if posslope_rel_phase[i] > np.pi:
        posslope_rel_phase[i] = posslope_rel_phase[i]%-np.pi
    elif posslope_rel_phase[i] < -1* np.pi:
        posslope_rel_phase[i] =  posslope_rel_phase[i]%np.pi
    else: continue

negslope_phase = avg_negslope_phase*np.pi/12               
negslope_meanphase = (np.arctan2(np.sin(negslope_phase).sum(),np.cos(negslope_phase).sum()))%(2*np.pi)
negslope_rel_phase = negslope_phase - negslope_meanphase

for i in range(len(negslope_rel_phase)):
    if negslope_rel_phase[i] > np.pi:
        negslope_rel_phase[i] = negslope_rel_phase[i]%-np.pi
    elif negslope_rel_phase[i] < -1* np.pi:
        negslope_rel_phase[i] =  negslope_rel_phase[i]%np.pi
    else: continue






xov_phase = round(pd.concat([posslope_phase, posslope_rel_phase, negslope_phase, negslope_rel_phase, peak_phase,
                       peak_rel_phase, trough_phase, trough_rel_phase], axis = 1)*(12/np.pi),2)
xov_phase.columns = ['XOV + Slope Pha (h)', 'XOV + Slope rPha (h)', 'XOV - Slope Pha (h)', 'XOV - Slope rPha (h)', 'XOV Peak Pha (h)', 
                     'XOV Peak rPha (h)', 'XOV Trough Pha (h)', 'XOV Trough rPha (h)']
xov_period = pd.concat([xov_posslope_period, xov_negslope_period, peaks_period , troughs_period], axis = 1 ).T
#groups by index
#xov_cycwise = xov_period.groupby(xov_period.index)
#calculates peak-trough apmlitude
xov_amp_cycwise = round(abs(peaks_val - troughs_val),4)
#mean peak-trough apm across cycles
xov_amp_mean = xov_amp_cycwise.mean(axis = 1) 


#calculates cyclewise  per mean of from all 4 methods
#xov_cycwise_per = xov_cycwise.mean().T


#calculates per across cycles
#xov_period = round(xov_cycwise_per.mean(axis = 1),2)



