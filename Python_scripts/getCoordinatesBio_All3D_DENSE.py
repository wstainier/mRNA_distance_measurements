# This code will create a series of CSV files which are necessary for performing the distance measurements in the following step. 
# These CSV files extract the coordinates of the germ granule segmentations from the Results files from the previous step. 
# Note that you might have to change the exact voxel sizes throughout the code to fit the specifications of your own image. This code was written for a voxel size of 
# 42.5nm x 42.5nm x 150nm (x,y,z). 

#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import csv
import pandas as pd 
import itertools
#from pyqtgraph.Qt import QtCore, QtGui
#import pyqtgraph.opengl as gl

'''
Reads and stores the intensity matrix for one z-slice 
Cuts off first row and first two columns
'''
def GetSliceIntensities(path):
    with open(path, 'r') as csv_file:
        matrix = []
        csv_reader = csv.reader (csv_file)
        for row in csv_reader:
            matrix.append(row)
        I = np.array(matrix)
        I = I[1:,2:] #cuts off first row and first two columns
        I = [[float(y) for y in x] for x in I]
        I = np.array(I)
    return I

'''
Returns x,y,z coordinates of fluorescence in current z-slice given the matrix of 
intensity values and which Z slice it is 
Assumes a voxel size of 42.5 x 42.5 x 150 nm (x, y, z)
'''
def getCoordinates(I,slice):
    y,x = I.nonzero() #y-rows, x-cols
    zlen = (len(y))
    z = [slice]*zlen
    x = np.array(x, dtype = float);y = np.array(y, dtype = float);z = np.array(z, dtype = float)
    #Scaling to biological size
    x = (x*42.5) + 21.25; y = (y*42.5) + 21.25; z = (z*150) + 75
    return(x,y,z)

'''
Similar to getCoordinates, but this function takes every voxel and divides it into 27 equal cuboids
The data points that are then extracted are the center points of every smaller cuboid
Assumes a voxel size of 42.5 x 42.5 x 150 nm (x, y, z)
'''

#Link for reference for np.meshgrid: https://stackoverflow.com/questions/1208118/using-numpy-to-build-an-array-of-all-combinations-of-two-arrays
def getCoordinatesDENSE(I,slice):
    y,x = I.nonzero() #y-rows, x-cols
    zlen = (len(y))
    z = [slice]*zlen
    x = np.array(x, dtype = float);y = np.array(y, dtype = float);z = np.array(z, dtype = float)
    #Scaling to biological size
    x = (x*42.5); y = (y*42.5); z = (z*150) 
    w = 42.5/(3*2); v = 150/(3*2)
    #Making coord the same dimensions as the slice coordinates to allow for concatenation 
    coord = np.zeros((1,3)) 
    coord = coord.astype(float)
    for i, j, k in zip(x, y, z):
            xlist = [i + w, i + 3*w, i + 5*w]
            ylist = [j + w, j + 3*w, j + 5*w]
            zlist = [k + v, k + 3*v, k + 5*v]
            mesh_store = np.array(np.meshgrid(xlist, ylist, zlist)).T.reshape(-1,3)
            coord = np.vstack((coord, mesh_store))
    coord = np.delete(coord,0,0)
    return(x,y,z, coord)

# Link for reference: https://stackoverflow.com/questions/64313602/numpy-opposite-of-nonzero-get-indices-of-zero-elements
#This function is to get all the zero points in the results CSV files --> these are the background
def getBackground(I, slice):
    y, x = np.nonzero(I == 0) #y-rows, x-cols
    zlen = (len(y))
    z = [slice]*zlen
    x = np.array(x, dtype = float);y = np.array(y, dtype = float);z = np.array(z, dtype = float)
    #Scaling to biological size
    x = (x*42.5) + 21.25; y = (y*42.5) + 21.25; z = (z*150) + 75
    return(x,y,z)


# In[2]:


#Defining the file (Note that you are already in the GranuleTranslation-Main folder)
userinput = input('Enter the name of the folder where the ZResults are stored:')
slicenumber = int(input('Number of slices in image '))


# In[16]:


#THIS IS FOR THE FILLED IN GRANULES
#Getting the coordinates for all the z-slices and storing it in an array called pos
pos = np.zeros((1,3)) #Making pos the same dimensions as the slice coordinates to allow for concatenation 
pos = pos.astype(float)
pos_back = np.zeros((1,3)) #This is for the background - i.e. the zero values 
pos_back = pos.astype(float)
for i in range(0,slicenumber,1):
    I = GetSliceIntensities(userinput + "/ZResults/Results"+str(i)+".csv")
    AxisLim = I.shape[0]#Number of rows/columns
    
    x,y,z = getCoordinates(I,i)
    SlicePos = np.dstack((x,y,z))
    SlicePos = SlicePos[0]
    SlicePos = SlicePos.astype(float)
    pos = np.vstack((pos, SlicePos))
    
    xback, yback, zback = getBackground(I,i)
    SlicePos_back = np.dstack((xback,yback,zback))
    SlicePos_back = SlicePos_back[0]
    SlicePos_back = SlicePos_back.astype(float)
    pos_back = np.vstack((pos_back, SlicePos_back))
pos = np.delete(pos,0,0) #Deleting the [0,0,0] used for initialization
pos_back = np.delete(pos_back,0,0) #doing the same for background 

#Checking to make sure that the background and the granule points are not overlapping using some size calculations
if (pos_back.size + pos.size)/(I.size*slicenumber) == 3:
    print('Size check good')
else: 
    print('ERROR! (in background and granule)')


# In[17]:


#Saving coordinates to file 
data = np.hsplit(pos,3)
X = np.array(data[0], dtype = float);Y = np.array(data[1], dtype = float);Z = np.array(data[2], dtype = float)
np.savetxt(userinput + '/3DCoordinates_Granule.csv', np.column_stack((X, Y, Z)), delimiter=",", fmt='%s')

data_back = np.hsplit(pos_back,3)
X_back = np.array(data_back[0], dtype = float);Y_back = np.array(data_back[1], dtype = float);Z_back = np.array(data_back[2], dtype = float)
np.savetxt(userinput + '/3DCoordinates_Background.csv', np.column_stack((X_back, Y_back, Z_back)), delimiter=",", fmt='%s')


# In[18]:

#THIS IS FOR THE 3D Surface (excluding edge cases)
#Getting the coordinates for all the z-slices and storing it in an array called pos_out
pos_out = np.zeros((1,3)) #Making pos_out the same dimensions as the slice coordinates to allow for concatenation 
pos_out = pos_out.astype(float)
pos_out_dense = np.zeros((1,3)) #same as above but for the "dense" arrays
pos_out_dense = pos_out_dense.astype(float)
for i in range(0,slicenumber,1):
    I = GetSliceIntensities(userinput + "/3DSurface_Results/Results"+str(i)+".csv")
    AxisLim = I.shape[0]#Number of rows/columns

    x,y,z = getCoordinates(I,i)
    size = len(z)
    SlicePos = np.dstack((x,y,z))
    SlicePos = SlicePos[0]
    SlicePos = SlicePos.astype(float)
    pos_out = np.vstack((pos_out, SlicePos))

    x_dense, y_dense, z_dense,coord_array = getCoordinatesDENSE(I,i)
    pos_out_dense = np.vstack((pos_out_dense, coord_array))

pos_out = np.delete(pos_out,0,0) #Deleting the [0,0,0] used for initialization
pos_out_dense = np.delete(pos_out_dense,0,0) #Deleting the [0,0,0] used for initialization

#Checking to make sure that the background and the granule points are not overlapping using some size calculations
if len(pos_out)*27 == len(pos_out_dense):
    print('Size check good – there are 27 times more "dense" points')
else: 
    print('ERROR! (in 3D Surface and 3D Surface DENSE)')

#Saving coordinates to file 
data = np.hsplit(pos_out,3)
X = np.array(data[0], dtype = float);Y = np.array(data[1], dtype = float);Z = np.array(data[2], dtype = float)
np.savetxt(userinput + '/3DCoordinates_3DSurface.csv', np.column_stack((X, Y, Z)), delimiter=",", fmt='%s')

data_dense = np.hsplit(pos_out_dense,3)
X_dense = np.array(data_dense[0], dtype = float);Y_dense = np.array(data_dense[1], dtype = float);Z_dense = np.array(data_dense[2], dtype = float)
np.savetxt(userinput + '/3DCoordinates_3DSurfaceDENSE.csv', np.column_stack((X_dense, Y_dense, Z_dense)), delimiter=",", fmt='%s')

#THIS IS FOR THE 3D Surface (including edge cases)
#Getting the coordinates for all the z-slices and storing it in an array called pos_out
pos_out = np.zeros((1,3)) #Making pos_out the same dimensions as the slice coordinates to allow for concatenation 
pos_out = pos_out.astype(float)
for i in range(0,slicenumber,1):
    I = GetSliceIntensities(userinput + "/3DSurface_Results_All/Results"+str(i)+".csv")
    AxisLim = I.shape[0]#Number of rows/columns
    x,y,z = getCoordinates(I,i)
    size = len(z)
    SlicePos = np.dstack((x,y,z))
    SlicePos = SlicePos[0]
    SlicePos = SlicePos.astype(float)
    pos_out = np.vstack((pos_out, SlicePos))
pos_out = np.delete(pos_out,0,0) #Deleting the [0,0,0] used for initialization

#Saving coordinates to file 
data = np.hsplit(pos_out,3)
X = np.array(data[0], dtype = float);Y = np.array(data[1], dtype = float);Z = np.array(data[2], dtype = float)
np.savetxt(userinput + '/3DCoordinates_3DSurface_All.csv', np.column_stack((X, Y, Z)), delimiter=",", fmt='%s')


# In[19]:


#THIS IS FOR THE Subtract ERODE --> it is from subtracting two images in ImageJ rather than running the Erode function
#Getting the coordinates for all the z-slices and storing it in an array called pos_erode
pos_erode = np.zeros((1,3)) #Making pos_erode the same dimensions as the slice coordinates to allow for concatenation 
pos_erode = pos_erode.astype(float)
for i in range(0,slicenumber,1):
    I = GetSliceIntensities(userinput + "/ErodeSubtract_Results/Results"+str(i)+".csv")
    AxisLim = I.shape[0]#Number of rows/columns
    x,y,z = getCoordinates(I,i)
    size = len(z)
    SlicePos = np.dstack((x,y,z))
    SlicePos = SlicePos[0]
    SlicePos = SlicePos.astype(float)
    pos_erode = np.vstack((pos_erode, SlicePos))
pos_erode = np.delete(pos_erode,0,0) #Deleting the [0,0,0] used for initialization

#Saving coordinates to file 
data = np.hsplit(pos_erode,3)
X = np.array(data[0], dtype = float);Y = np.array(data[1], dtype = float);Z = np.array(data[2], dtype = float)
np.savetxt(userinput + '/3DCoordinates_ErodeSubtract.csv', np.column_stack((X, Y, Z)), delimiter=",", fmt='%s')

#Checking that the size of all granule points is equal to the outline points + the erode points 
if pos.size == pos_out.size + pos_erode.size:
    print('Erode and outline equal to granules')
else: 
    print('ERROR! (in erode and outline)')
    print("Granule: " + str(pos.size) + "    3D surface: " + str(pos_out.size) + "    erode subtract: " + str(pos_erode.size) + ".")




###This next bit of code will get the centroid data in the same format as the above data (CSV files with 3 columns – position in biologically scaled
### size)

#For the centroid data that EXCLUDES the edge cases

df = pd.read_csv(userinput + "/CentroidData/3D_ObjectCount_ResultsTable.csv", usecols=["X", "Y", "Z"])
df["X"] = df["X"] * 42.5 + 21.25
df["Y"] = df["Y"] * 42.5 + 21.25
df["Z"] = df["Z"] * 150 + 75

centroid_size = df.size

df.to_csv(userinput + '/Centroid.csv', header = False, index = False)

#For the centroid data that INCLUDES the edge cases

df = pd.read_csv(userinput + "/CentroidData_All/3D_ObjectCount_ResultsTable.csv", usecols=["X", "Y", "Z"])
df["X"] = df["X"] * 42.5 + 21.25
df["Y"] = df["Y"] * 42.5 + 21.25
df["Z"] = df["Z"] * 150 + 75

centroidall_size = df.size

df.to_csv(userinput + '/Centroid_All.csv', header = False, index = False)

#Checking that the centroid dataframes of excluding and including the edge cases are different 
if centroid_size == centroidall_size:
    print('ERROR! The INCLUDE and EXCLUDE edge cases centroid dataframes are of equal size')
else: 
    print('The two different centroid dataframes (include and exclude edge cases) are different, as expected.')





