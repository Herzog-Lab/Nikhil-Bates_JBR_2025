# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:04:15 2021

@author: Nikhil
"""

startTime = datetime.now()

filename_raw = input("Enter filename without extention(.tif):")
filename_full = str(filename_raw) + str(file_extention)

sampling_interval_mins = input("Enter sampling interval (time in MINS between consecutive samples): ")
sampling_interval_mins = int(sampling_interval_mins)

#creates folders to save results
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

createFolder('./Results_' +str(filename_raw))
createFolder('./Results_' +str(filename_raw)+'/cache')
createFolder('./Results_' +str(filename_raw)+'/Plots')
createFolder('./Results_' +str(filename_raw)+'/Plots/Maps')
createFolder('./Results_' +str(filename_raw)+'/Plots/Maps/Period Maps')
createFolder('./Results_' +str(filename_raw)+'/Plots/Maps/Phase Maps')
createFolder('./Results_' +str(filename_raw)+'/Plots/Maps/Amp Maps')
createFolder('./Results_' +str(filename_raw)+'/Plots/Stack Plots')
createFolder('./Results_' +str(filename_raw)+'/Plots/Stack Plots/COM')
createFolder('./Results_' +str(filename_raw)+'/Plots/Stack Plots/CWT Abs Phase Evol Map')
createFolder('./Results_' +str(filename_raw)+'/Plots/Stack Plots/CWT Rel Phase Evol Polar')
createFolder('./Results_' +str(filename_raw)+'/Plots/Stack Plots/CWT Abs Phase Evol Polar')
createFolder('./Results_' +str(filename_raw)+'/Plots/Stack Plots/Polar Histograms')
createFolder('./Results_' +str(filename_raw)+'/Plots/Stack Plots/Polar Histograms/CWT Abs Phase Polar Histograms')

#set bin size for histograms
bin_size = 10

#reads image
image = io.imread(filename_full)

#measure pixel intensities
pixs = []
for i in range(image.shape[1]):
    for j in range(image.shape[2]):
        pixs.append(image[:, i,j])

#sums values along time axis
pixs_sum = np.sum(np.array(pixs),axis = 1)
#filters out pixels with total signal = 0 and counts total pixels with signal
valid_pixs = len(pixs_sum[pixs_sum >0])


ratio = round(pixels_per_image/valid_pixs,2)
 #sets scaling factor based on total number of pixels with signal. If time taken for analysis is too large, reduce the pixels_per_image values

if ratio <= 0.05:
    image_scaling_factor = 0.1
elif 0.05 < ratio <= 0.1:
    image_scaling_factor = 0.2
elif 0.1 < ratio <= 0.16:
    image_scaling_factor = 0.3        
elif 0.16 < ratio <= 0.25:
    image_scaling_factor = 0.4
elif 0.25 < ratio <= 0.37:
    image_scaling_factor = 0.5
elif 0.37 < ratio <= 0.49:
    image_scaling_factor = 0.6
elif 0.49 < ratio <= 0.64:
    image_scaling_factor = 0.7        
elif 0.64 < ratio <= 0.82:
    image_scaling_factor = 0.8
else:
    image_scaling_factor = 0.0

if smooth == 'default':
        smooth = 1
else: 
    smooth = smooth




#rescales image based on above conditions
dim = rescale(image, image_scaling_factor)
im = np.empty((image.shape[0],dim.shape[1],dim.shape[2]))
for i in range(image.shape[0]):
    im[i] = rescale(image[i], image_scaling_factor)
im = gaussian(im, sigma=smooth)


io.imshow(im[0])
plt.show()
plt.close()