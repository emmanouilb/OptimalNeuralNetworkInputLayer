from PIL import Image

from scipy import misc
import numpy as np


im = Image.open("testing_set_Parlitz_KS.ppm")
im.show()

#

from scipy.misc import imread
img_testing_set = imread("testing_set_Parlitz_KS.ppm")
img_testing_set.shape
img_training_set = imread("training_set_Parlitz_KS.ppm")
img_training_set.shape
img_forecast_set = imread("forecast_set_Parlitz_KS.ppm")
img_forecast_set.shape

import matplotlib.pyplot as plt # import
plt.imshow(img_training_set) #load
plt.show()  # show the window

def average(pixel):
    return (pixel[0] + pixel[1] + pixel[2]) / 3
    

import numpy as np
img_training_set_grey = np.zeros((img_training_set.shape[0], img_training_set.shape[1])) # init 2D numpy array
# get row number
for rownum in range(len(img_training_set)):
    for colnum in range(len(img_training_set[rownum])):
        img_training_set_grey[rownum][colnum] = average(img_training_set[rownum][colnum])

import matplotlib.cm as cm # 

plt.imshow(img_training_set_grey, cmap = cm.Greys_r)
plt.show()

# problems, it does not seem like the right picture.
#misc.imsave('img_training_set_grey.ppm', img_training_set_grey)

# try again after using i_view32.exe for greyscale

img_testing_set_grey = imread("testing_set_Parlitz_KS_grey.ppm")
img_testing_set_grey.shape
img_training_set_grey = imread("training_set_Parlitz_KS_grey.ppm")
img_training_set_grey.shape
img_forecast_set_grey = imread("forecast_set_Parlitz_KS_grey.ppm")
img_forecast_set.shape

# show they are all the same now
np.average(img_training_set_grey[:,:,0])
np.average(img_training_set_grey[:,:,1])
np.average(img_training_set_grey[:,:,2])

# extract Fig3a from 
# /cygdrive/c/Users/eurico/SunspotAnalysis/Papers/PhysRevLett.84.1890.pdf
data_img_testing_set_grey=img_testing_set_grey[:,:,0]
data_img_training_set_grey=img_training_set_grey[:,:,0]
data_img_forecast_set_grey=img_forecast_set_grey[:,:,0]
data_img_testing_set_grey=data_img_testing_set_grey/255.0*1.5 # according to paper
data_img_training_set_grey=data_img_training_set_grey/255.0*1.5 # according to paper
data_img_forecast_set_grey=data_img_forecast_set_grey/255.0*1.5 # according to paper
from scipy import stats
data_img_testing_set_grey.shape
data_img_training_set_grey.shape
data_img_forecast_set_grey.shape
stats.describe(data_img_testing_set_grey.flatten())
stats.describe(data_img_training_set_grey.flatten())
stats.describe(data_img_forecast_set_grey.flatten())

data_img_training_set_grey_histogram = np.histogram(data_img_training_set_grey.flatten(), bins=np.arange(0,1.5,0.1), density=True)
plt.hist(data_img_training_set_grey.flatten(), bins='auto')
plt.show()

# export to csv
np.savetxt("data_img_testing_set_grey.csv", data_img_testing_set_grey, delimiter=",")
np.savetxt("data_img_training_set_grey.csv", data_img_training_set_grey, delimiter=",")
np.savetxt("data_img_forecast_set_grey.csv", data_img_forecast_set_grey, delimiter=",")












