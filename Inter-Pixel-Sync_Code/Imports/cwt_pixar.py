# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 18:07:26 2021

@author: Nikhil
"""

#CWT analysis
#CWT analysis works better with smooth waveforms. So we will use a higher sigma gaussian filter
#to smooth pixel values
#reads image 
image = io.imread(filename_full)
dim = rescale(image, image_scaling_factor, anti_aliasing=True)
im = np.empty((image.shape[0],dim.shape[1],dim.shape[2]))
for i in range(image.shape[0]):
    im [i] = rescale(image[i], image_scaling_factor, anti_aliasing=True)
im = gaussian(im, sigma=smooth) #- gaussian(im, sigma=20.0)

#reads gaussian filtered (smoothed) values for pixels reported as rhythmic by metacycle
pix_val_ts = []
for i in range(len(rhythmic_pixs)):
    pix_val_ts.append(im[:,rhythmic_pixs['Y'].loc[i],rhythmic_pixs['X'].loc[i]])
pix_val_ts = pd.DataFrame(pix_val_ts).T/pd.DataFrame(pix_val_ts).mean(axis = 1)

periods = np.round(np.linspace(min_per*60, (max_per)*60, nper), 2)
cwt_period = []
cwt_phase = []
cwt_power = []
detrended_ts = []
normalized_ts = []
FFT_spectrum = []
amp_envelope = []
wAn = WAnalyzer(periods, sampling_interval_mins, time_unit_label = "min")
for i in range(len(pix_val_ts.columns)):
    filtered = wAn.sinc_detrend(signal = pix_val_ts[i], T_c = T_c*60)
    normalized = wAn.normalize_amplitude(filtered, 24*60)
    modulus, transform = wAn.compute_spectrum(normalized, do_plot=False)
    FFT = wAn.get_averaged_spectrum()
    ridge_results = wAn.get_maxRidge()
    amp =  wAn.get_envelope(filtered, window_size = 24*60, SGsmooth = True)
    cwt_period.append(ridge_results["periods"]/60)
    cwt_phase.append(ridge_results["phase"])
    cwt_power.append(ridge_results["power"])
    detrended_ts.append(filtered)
    normalized_ts.append(normalized)
    FFT_spectrum.append(FFT)
    amp_envelope.append(amp)

FFT_maxpower_index = np.argmax(FFT_spectrum, axis = 1)
FFT_period = []
for i in FFT_maxpower_index:
    FFT_period.append(round(periods[i]/60, 2))
FFT_period = pd.DataFrame(FFT_period)

#smoothing instantaneous period
pad = int((60/sampling_interval_mins)*6)
roll_win = int(pad*2)
cwt_per= pd.DataFrame(cwt_period).T
cwt_fpad = pd.concat([pd.DataFrame(cwt_per[:pad].mean()).T]*pad)
cwt_bpad = pd.concat([pd.DataFrame(cwt_per[len(cwt_per)-pad:].mean()).T]*pad)
cwt_per_pad = pd.concat([cwt_fpad, cwt_per, cwt_bpad])
cwt_period = cwt_per_pad.rolling(roll_win, center = True).mean()
cwt_period = cwt_period[pad:len(cwt_period)-pad]
avg_cwt_per = pd.DataFrame(round(cwt_period[pad:len(cwt_period)-pad].mean(),2))
avg_inst_per = round(cwt_period.mean(axis = 1),2)
#avg_inst_per_std =round(cwt_period.std(),2)

avg_cwt_amp = pd.DataFrame(amp_envelope).T
avg_cwt_amp = avg_cwt_amp.loc[(60/sampling_interval_mins)*6 : len(amp_envelope) -
                (60/sampling_interval_mins)*6].mean().reset_index(drop = True)

avg_ts = pd.DataFrame(detrended_ts).T.mean(axis = 1)
avg_inst_amp = pd.DataFrame(amp_envelope).T.mean(axis =1)

cwt_phase = pd.DataFrame(cwt_phase).reset_index(drop = True)
avg_inst_phase = np.arctan2(np.sin(cwt_phase).sum(axis=0), np.cos(cwt_phase).sum(axis=0))%(2*np.pi)
avg_inst_phase_h = round(pd.DataFrame(avg_inst_phase *12/np.pi),2) 
cwt_phase_h = round(pd.DataFrame(cwt_phase) * 12 / np.pi,2)


#Kuromoto Order paramater analysis
#first order kuromoto order parameter
r1_df1 = cwt_phase.dropna(axis = 0)
r1 = pd.DataFrame()
for column in cwt_phase:
    r1 = r1.append([sum((np.e ** (1j * r1_df1[column]))/len(r1_df1))])
r1 = abs(r1).reset_index(drop = True)
r1.columns = ['First Order SI']

# second order kuromoto order parameter
r2 = pd.DataFrame()
r2_phases = 2*r1_df1
for column in r2_phases:
    r2 = r2.append([sum((np.e ** (1j * r2_phases[column]))/len(r2_phases))])
r2 = abs(r2).reset_index(drop = True)
r2.columns = ['Sec Order SI']






#saves CWT data
cwt_period.to_csv("./Results_" +str(filename_raw)+"/CWT Instantaneous Period.csv", 
                                  header= np.arange(1, len(cwt_period.columns)+1))
pd.DataFrame(cwt_phase).T.to_csv("./Results_" +str(filename_raw)+"/CWT Instantaneous Phase_rad.csv", 
                                  header= np.arange(1, len(cwt_phase)+1))
cwt_phase_h.T.to_csv("./Results_" +str(filename_raw)+"/CWT Instantaneous Phase_h.csv", 
                                  header= np.arange(1, len(cwt_phase)+1))
pd.DataFrame(cwt_power).T.to_csv("./Results_" +str(filename_raw)+"/CWT Instantaneous Period Power.csv", 
                                  header= np.arange(1, len(cwt_power)+1))
pd.DataFrame(detrended_ts).T.to_csv("./Results_" +str(filename_raw)+"/Detrended Traces.csv", 
                                  header= np.arange(1, len(detrended_ts)+1))
pd.DataFrame(normalized_ts).T.to_csv("./Results_" +str(filename_raw)+"/Normalized Traces.csv", 
                                  header= np.arange(1, len(normalized_ts)+1))
pd.DataFrame(FFT_spectrum).T.to_csv("./Results_" +str(filename_raw)+"/FFT Spectrum.csv", 
                                  header= np.arange(1, len(FFT_spectrum)+1))
pd.DataFrame(amp_envelope).T.to_csv("./Results_" +str(filename_raw)+"/CWT Instantaneous Amplitude.csv", 
                                  header= np.arange(1, len(amp_envelope)+1))

FFT_period.to_csv("./Results_" +str(filename_raw)+"/FFT Period.csv", header = False)

ensemble_dyn = pd.concat([avg_inst_per, avg_inst_amp, avg_inst_phase_h, round(r1,2), round(r2,2)], axis = 1)
ensemble_dyn.columns = ['Per (h)', 'Amp', 'Pha (h)', 'Sync Index 1', 'Sync Index 2']
ensemble_dyn.to_csv("./Results_" +str(filename_raw)+"/Ensemble Dynamics.csv")
