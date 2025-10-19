


#measure intensities and creates list of  roi coordinates and time points
x = []
y = []
pix_val_ts = []
for i in range(im.shape[1]):
    for j in range(im.shape[2]):
        pix_val_ts.append(im[:, i,j])
        y.append(i)
        x.append(j)

xy = pd.concat([pd.DataFrame(x), pd.DataFrame(y)], axis = 1)
xy.columns = ['X', 'Y']
CycID = pd.DataFrame(np.arange(1, len(xy)+1), columns = ['CycID'])
time = [np.arange(0, im.shape[0]*sampling_interval_mins, sampling_interval_mins)]

ref_time = pd.DataFrame(time).T
xticks = np.round(ref_time/60,1)
ncycles = int(len(ref_time[0])*sampling_interval_mins/(24*60))
vlines= np.arange(0, int((24*ncycles*60/sampling_interval_mins)+1), int(6*60/sampling_interval_mins))/(60/sampling_interval_mins)

#plots  time series for randomly selected raw timeseries
if len(pix_val_ts) > 100:
    n =100
else:
    n = int(len(pix_val_ts)/2)

for i in range(n):
    sample = pix_val_ts[random.randint(0, len(pix_val_ts)-1)]
    if sample.mean() != 0:
        plt.plot(xticks[0], sample/sample.mean())
    else:
        continue

for k in vlines:
    plt.axvline(x=k, color='k', linestyle='--', alpha = 0.2)
plt.title("Representative Raw Time Series (sigma = " + str(smooth) + ")")
plt.xlabel("Time (h)")
plt.ylabel("Pixel Intensity (mean normalized)")
plt.savefig("./Results_" +str(filename_raw)+ "/Plots/Sample Time Series.tif", format="TIF", dpi=300)
plt.show()
plt.close()


phasemap_time = float(input("Enter time (IN HOURS) for phasemap generation : "))
while phasemap_time*60  not in time[0]*1.0:
    warnings.warn("Time entered for phasemap generation is either out of range or not compatible with sampling interval. Please enter a different time")
    phasemap_time = float(input("Enter time (IN HOURS) for phasemap generation : "))

print("COM analysis in progress..")
#COM analysis
image_dimension = pd.DataFrame([image.shape[2], image.shape[1]]).T
image_dimension.columns = ['X', 'Y']
image_dimension.to_csv('./Results_'+str(filename_raw)+'/COM Image Dimension.csv')
com_x = []
com_y = []
for i in range(image.shape[0]):
    com = ndimage.measurements.center_of_mass(image[i])
    com_x.append(com[0])
    com_y.append(com[1])
com = pd.DataFrame(list(zip(com_x,com_y)), columns = ['com_x', 'com_y'])
com.to_csv("./Results_" +str(filename_raw)+"/COM cooridnates.csv")


print("COM analysis completed. Metacycle Analysis in progress..")
#Metacycle Analysis
#creates and saves pixelwise timeseries csv file as input metacycle analysis
metadata= pd.DataFrame(pix_val_ts).reset_index(drop = True)
metadata_colsum = pd.DataFrame(metadata.sum(axis = 1), columns = ['colsum'])
metadata_filt = pd.concat([xy, CycID, metadata, metadata_colsum], axis = 1)
metadata_filt= metadata_filt[metadata_filt['colsum'] > 0].reset_index(drop=True)
xy = metadata_filt.iloc[:,0:2]
metadata_norm= metadata_filt.iloc[:, 3:len(metadata_filt.columns)-1].T
metadata_norm = metadata_norm/metadata_norm.mean().T
metadata_norm = metadata_norm.T
metadata_norm.columns =  time
metadata_norm.insert(0, "CycID", metadata_filt['CycID'], True)
metadata_norm.to_csv("metadata.csv", index = False)


#sources to R to run metacycle analysis
if platform.system() == 'Windows':
    r_source = ro.r['source']
    r_source('./Imports/meta2d_script.R')
elif platform.system() == 'Darwin':
    subprocess.run (['/usr/local/bin/Rscript', './Imports/meta2d_script.R'])

#reads and filters arrhythmic pixels data based on metacycle BH_Q index
rhythmic_pixs = pd.read_csv("./metaout/meta2d_metadata.csv")
rhythmic_pixs = pd.concat([xy,rhythmic_pixs], axis = 1)
rhythmic_pixs = rhythmic_pixs[rhythmic_pixs['meta2d_BH.Q'] < 0.01].reset_index(drop = True)
corrected_metaphase = pd.DataFrame(((rhythmic_pixs['meta2d_phase']-720))%(24*60))
rhythmic_pixs = pd.concat([rhythmic_pixs['X'], rhythmic_pixs['Y'], rhythmic_pixs['CycID'],
                           round(rhythmic_pixs['meta2d_period']/60, 2), 
                           round(corrected_metaphase/60,2),
                           rhythmic_pixs['meta2d_rAMP'],rhythmic_pixs['meta2d_BH.Q']], axis = 1)

#saves arrhythmic cells as a separate file
arrhythmic_pixs = pd.read_csv("./metaout/meta2d_metadata.csv")
arrhythmic_pixs = pd.concat([xy,arrhythmic_pixs], axis = 1)
arrhythmic_pixs = arrhythmic_pixs[arrhythmic_pixs['meta2d_BH.Q'] > 0.01].reset_index(drop = True)
corrected_ar_metaphase = pd.DataFrame(((arrhythmic_pixs['meta2d_phase']-720))%(24*60))
arrhythmic_pixs = pd.concat([arrhythmic_pixs['X'], arrhythmic_pixs['Y'], arrhythmic_pixs['CycID'],
                           round(arrhythmic_pixs['meta2d_period']/60, 2), arrhythmic_pixs['meta2d_BH.Q']], axis = 1)

#reads image 
image = io.imread(filename_full)
dim = rescale(image, image_scaling_factor, anti_aliasing=True)
im = np.empty((image.shape[0],dim.shape[1],dim.shape[2]))
for i in range(image.shape[0]):
    im [i] = rescale(image[i], image_scaling_factor, anti_aliasing=True)
im = gaussian(im, sigma=smooth) #- gaussian(im, sigma=20.0)

#reads gaussian filtered (smoothed) values for pixels reported as rhythmic by metacycle
ar_pix_val_ts = []
for i in range(len(arrhythmic_pixs)):
    ar_pix_val_ts.append(im[:,arrhythmic_pixs['Y'].loc[i],arrhythmic_pixs['X'].loc[i]])
ar_pix_val_ts = pd.concat([arrhythmic_pixs, pd.DataFrame(ar_pix_val_ts)], axis = 1)




percnt_rhy = pd.DataFrame([len(metadata_norm), len(rhythmic_pixs), 
                           round((len(rhythmic_pixs)/(len(metadata_norm)))*100, 2)], 
                          index = ["Total pixels with signal", "Total Rhythmic pixels", "Percent Rhythmic"])

percnt_rhy.to_csv("./Results_" +str(filename_raw)+"/Percent Rhythmic.csv", header = False)
