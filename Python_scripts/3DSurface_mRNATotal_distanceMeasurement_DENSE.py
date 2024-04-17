#This code will perform the analysis of the points which are categorized as total mRNA. It will output a distance measurement of those points to the surface 
#of the granule.
#If needed, you will need to change the value that is half the length of the diagonal of the voxel in the last for loop of the file. It is currently set at 81. 

#!/usr/bin/env python
# coding: utf-8

# In[185]:


#Find the sub-pixel location of total mRNA (from FISH-Quant)
#Find the closest 3D surface point to each mRNA point and plot distances
#Plot the center point, granules, and closest point on the granule

import numpy as np
import csv
import pandas as pd
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

def getPoints(filename):
    x = list(); y = list(); z = list()
    with open (filename, 'r', encoding='utf-8-sig') as csv_file:
        csv_reader = csv.reader (csv_file)
        for line in csv_reader:
                x.append(line[0]); y.append(line[1]); z.append(line[2])
    x = np.array(x); y = np.array(y); z = np.array(z)
    x = x.astype(float);y = y.astype(float);z = z.astype(float)
    arr = np.stack((x, y, z), axis=-1)
    arr = arr.astype(float)
    return (arr)


# In[186]:


# #Defining the file (Note that you are already in the GranuleTranslation-Main folder)
userinput = input('Enter the name of the folder where the 3DCoordinates files are stored:')    

arr1 = getPoints(userinput + '/3DCoordinates_3DSurfaceDENSE.csv') #3D surface excluding edges – DENSE points
arr2 = getPoints(userinput + '/3DCoordinates_3DSurface.csv')#3D surface – will leave this named as "closest_outline" etc. in the code 
arr3 = getPoints(userinput + '/FQ3DCoordinatesTotalmRNA.csv')#total mRNA
arr4 = getPoints(userinput + '/3DCoordinates_Granule.csv')#each pixel of granule
arr5 = getPoints(userinput + '/3DCoordinates_ErodeSubtract.csv') #the granule points which are NOT on edge – based on image subtraction
arr6 = getPoints(userinput + '/3DCoordinates_Background.csv') #the background points
arr7 = getPoints(userinput + '/3DCoordinates_3DSurface_All.csv') #3D surface ALL


# In[188]:


#STEP 2 - Find the outline point and its index that each mRNA total point is closest to. Also find the granule point and its index that each mRNA total point is closest to
closest_outline_dense = list() #List of the closest outline point – this is really the DENSE 3D surface 
minindex_dense = list() #List of closest outline index – this is really the DENSE 3D surface
closest_outline = list() #List of the closest outline point – this is really the 3D surface 
minindex = list() #List of closest outline index – this is really the 3D surface
closest_gra = list() #List of the closest granule point
minindex_gra = list() #List of closest granule index
closest_erode = list() #List of the closest eroded granule point
minindex_erode = list() #List of closest eroded granule index
closest_back = list() #List of the closest background point
minindex_back = list() #List of closest background index
closest_3D_All = list() #List of the closest 3D Surface ALL point
minindex_3D_All = list() #List of closest 3D Surface ALL point
for c in arr3:
    curindex_dense = cdist([c], arr1).argmin()
    minindex_dense.append(curindex_dense)
    closest_outline_dense.append(arr1[curindex_dense])

    curindex = cdist([c], arr2).argmin()
    minindex.append(curindex)
    closest_outline.append(arr2[curindex])
    
    curindex_gra = cdist([c], arr4).argmin()
    minindex_gra.append(curindex_gra)
    closest_gra.append(arr4[curindex_gra])
    
    curindex_erode = cdist([c], arr5).argmin()
    minindex_erode.append(curindex_erode)
    closest_erode.append(arr5[curindex_erode])
    
    curindex_back = cdist([c], arr6).argmin()
    minindex_back.append(curindex_back)
    closest_back.append(arr6[curindex_back])

    curindex_3D_All = cdist([c], arr7).argmin()
    minindex_3D_All.append(curindex_3D_All)
    closest_3D_All.append(arr7[curindex_3D_All])

minindex_dense = np.array(minindex_dense, dtype = int)
closest_outline_dense = np.array(closest_outline_dense, dtype = float)

minindex = np.array(minindex, dtype = int)
closest_outline = np.array(closest_outline, dtype = float)

minindex_gra = np.array(minindex_gra, dtype = int)
closest_gra = np.array(closest_gra, dtype = float)

minindex_erode = np.array(minindex_erode, dtype = int)
closest_erode = np.array(closest_erode, dtype = float)

minindex_back = np.array(minindex_back, dtype = int)
closest_back = np.array(closest_back, dtype = float)

minindex_3D_All = np.array(minindex_3D_All, dtype = int)
closest_3D_All = np.array(closest_3D_All, dtype = float)


#STEP 3 - Find and plot distance between closest 3D surface point and mRNA total points & categorize whether point is inside or outside granule
mindist_dense = list()
for c2,translation in zip(closest_outline_dense, arr3): 
    mindist_dense.append(np.linalg.norm(c2-translation))
mindist_dense = np.array(mindist_dense, dtype = float)
np.savetxt(userinput + '/mRNATotal_3DSurfaceDENSE_Distance.csv', mindist_dense ,delimiter=",", fmt='%s')
print('DENSE 3D Surface CSV saved')

mindist = list()
for c2,translation in zip(closest_outline, arr3): 
    mindist.append(np.linalg.norm(c2-translation))
mindist = np.array(mindist, dtype = float)
np.savetxt(userinput + '/mRNATotal_3DSurface_Distance.csv', mindist ,delimiter=",", fmt='%s')
print('3D Surface CSV saved')

#This is the list to the closest granule point
mindist_gra = list()
for c2,translation in zip(closest_gra, arr3): 
    mindist_gra.append(np.linalg.norm(c2-translation))
mindist_gra = np.array(mindist_gra, dtype = float)
np.savetxt(userinput + '/mRNATotal_Granule_Distance.csv', mindist_gra ,delimiter=",", fmt='%s')
print('Granule CSV saved')

#This is the list to the closest eroded granule point
mindist_erode = list()
for c2,translation in zip(closest_erode, arr3): 
    mindist_erode.append(np.linalg.norm(c2-translation))
mindist_erode = np.array(mindist_erode, dtype = float)
np.savetxt(userinput + '/mRNATotal_ErodeSubtract_Distance.csv', mindist_erode ,delimiter=",", fmt='%s')
print('Eroded (subtract) Granule CSV saved')

#This is the list to the background point
mindist_back = list()
for c2,translation in zip(closest_back, arr3): 
    mindist_back.append(np.linalg.norm(c2-translation))
mindist_back = np.array(mindist_back, dtype = float)
np.savetxt(userinput + '/mRNATotal_Background_Distance.csv', mindist_back ,delimiter=",", fmt='%s')
print('Background Granule CSV saved')

#This is the list to the closest 3D Surface All granule point
mindist_3D_All = list()
for c2,translation in zip(closest_3D_All, arr3): 
    mindist_3D_All.append(np.linalg.norm(c2-translation))
mindist_3D_All = np.array(mindist_3D_All, dtype = float)
np.savetxt(userinput + '/mRNATotal_3DSurfaceALL_Distance.csv', mindist_3D_All ,delimiter=",", fmt='%s')
print('3D Surface ALL Granule CSV saved')

##Finding the mRNA points that are closest to 3D surfaces which are in both datasets

df_trans = pd.read_csv(userinput + '/FQ3DCoordinatesTotalmRNA.csv', header = None, names = ['X', 'Y', 'Z'])

df_all = pd.read_csv(userinput + '/mRNATotal_3DSurfaceALL_Distance.csv', header = None, names = ['All_3D_Surface'])

df_3D = pd.read_csv(userinput + '/mRNATotal_3DSurface_Distance.csv', header = None, names = ['3D_Surface'])

df_3D_DENSE = pd.read_csv(userinput + '/mRNATotal_3DSurfaceDENSE_Distance.csv', header = None, names = ['DENSE_3D_Surface'])

both_df = [df_all, df_3D]

df = pd.concat(both_df, axis = 1)

#Calculating the difference between the two distance values. If the distance values are equal, the difference should
#be 0. These are thus the points that were closest to a 3D surface that was in both the "Edge excluded" and the 
#"Edge included" dataset. 
df['Difference'] = df['All_3D_Surface'] - df['3D_Surface'] 

#Quick check to make sure that the minimum distance to the "edge excluded" 3D Surface is always smaller or equal to  
#the "edge included" 3D Surface.
if sum(df['Difference'] > 0) == 0:
    print('As expected, all the values of the difference are smaller or equal to 0.')
else: 
    print('ERROR! There are some difference values that are greater than 0.')

consistent_index = df.index[df['Difference'] == 0].tolist() #finding the index of where the difference is equal to 0

#getting the location of mRNA total points whose closest DENSE 3D Surface is in both the "edge included" and "edge excluded"
#datasets – also concatenating the distance from DENSE 3D Surface to mRNA total point 
new_trans_df = df_trans.iloc[consistent_index]  
new_trans_df_dist = df_3D_DENSE.iloc[consistent_index]  

both_trans_df = [new_trans_df, new_trans_df_dist]

final_trans_df = pd.concat(both_trans_df, axis = 1)

final_trans_df = final_trans_df.rename(columns = {'3D_Surface':'Distance_to_DENSE_3D_Surface'})

#Now use the .to_csv to change that back to csv to get the final data – this CSV stores the x,y,z coordinates of the mRNA total points used 
#as well as the distance to the DENSE 3D surface (distance is 4th column)

final_trans_df.to_csv(userinput + '/3D_Surface_Consistent_and_mRNATotalDistance.csv', header = True, index = False)
print('The number of consistent mRNA total points is '  + str(len(consistent_index)) + ' out of ' + str(len(df_trans)) + ' total points.')


###Finding the location of the points (inside vs outside the granule) - but adjusting the distance to the DENSE 3D surface to either negative or positive

location_list = np.empty(len(mindist), dtype=np.dtype('U100'))
adjusted_mindist = np.zeros(len(mindist))
for i in range(0,len(mindist)):
    if mindist_gra[i] < 81: #really it is the half length of the diagonal of the cuboid that is pixel – can put exact number
        if mindist_gra[i] == mindist[i]:
            if mindist_erode[i] < mindist_back[i]:
                location_list[i] = 'inside edge'
                adjusted_mindist[i] = -1*mindist_dense[i]
            elif mindist_erode[i] > mindist_back[i]:
                location_list[i] = 'outside edge'
                adjusted_mindist[i] = mindist_dense[i]
            elif mindist_erode[i] == mindist[i]:
                print('Error at ' + str([i]) + '? The subtract erode and 3D surface minimum distance are the same.')
            else: 
                print('Error at ' + str([i]) + '? Maybe negative values?')
                print('Mindist_erode' + str([i]) + str(mindist_erode[i]))
                print('Mindist_back' + str([i]) + str(mindist_back[i]))
        else: 
            location_list[i] = 'inside'
            adjusted_mindist[i] = -1*mindist_dense[i]
    else:
        location_list[i] = 'outside'
        adjusted_mindist[i] = mindist_dense[i]


df = pd.DataFrame({"Distance to DENSE 3D Surface": mindist_dense, "Location": location_list, "Adjusted Distance to DENSE 3D Surface": adjusted_mindist})
df.to_csv(userinput + '/FullLocation_mRNATotal_Outline_Distance_notConsistent.csv', index=False)

#This will be for all points (not just the consistent ones?)
unique, counts = np.unique(location_list, return_counts=True)
location_dict = dict(zip(unique, counts))
print('The percentage of edge cases is: ' + str((location_dict['inside edge'] + location_dict['outside edge'])/len(location_list)*100) + '.')

df_consistent_adjusteddist = df.iloc[consistent_index]

all_trans_df = [final_trans_df, df_consistent_adjusteddist]

final_all_trans_df = pd.concat(all_trans_df, axis = 1)

#Now use the .to_csv to change that back to csv to get the final data – this CSV stores the x,y,z coordinates of the mRNA total points used 
#as well as the distance to the 3D surface (distance is 4th column)

final_all_trans_df.to_csv(userinput + '/ALL_DENSE_3D_Surface_Consistent_and_mRNATotalDistance.csv', header = True, index = False)


















