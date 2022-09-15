"""More patients
"""
import os
import numpy as np
import pydicom
from matplotlib import pyplot as plt
import cv2
import funcs as fc
import pandas as pd
import paths as pa #import paths
from skimage.draw import polygon

#import contours_1 as cnt1
#import inside_pol as ip
from shapely import geometry

import min_distance_funcs as mdf

"""Defining files and folders path"""

ids_list=pa.ids_list

patients=list(pd.read_csv(ids_list)["PatientId"])

plots_folder=pa.plots_folder

min_dist_ind_pats_folder=pa.min_dist_ind_pats_folder

min_dist_df=pa.min_dist_df

#pa.main_path


"""Defining functions"""

def get_Body_ROI_names(RS_path):

	rs = pydicom.read_file(RS_path)

	body_names=[]
	for key in rs.StructureSetROISequence:
		
		ROI=key.ROIName

		#print(ROI)

		if (ROI[0:4]=="BODY") | (ROI[0:4]=="Body"):

			body_names.append(ROI)

	return list(np.unique(body_names))




def min_dist_all_pat(patients):



	df_row=pd.DataFrame(columns = ['z_inter','min_dist','ID'])
	for pat_id in patients:

		path=pa.main_path+str(pat_id)+'/'

		RS_path=path+fc.find_RS_files(path)

		print("\n\n PATIENT ID",pat_id)

		body_names=get_Body_ROI_names(RS_path)

		print(body_names)

		print("getting PTV contours")

		x_PTV,y_PTV,z_PTV=fc.get_contour_all_slice("PTV_ALL",RS_path)

		x_PTV,y_PTV,z_PTV=mdf.sort_contours(x_PTV,y_PTV,z_PTV)


		df=mdf.get_cont_all_fractions(body_names,RS_path,x_PTV,y_PTV,z_PTV,str(pat_id))

		df.to_csv(min_dist_ind_pats_folder+"min_dist_pat"+str(pat_id))

		df_row = pd.concat([df_row, df])

	return df_row

"""Runnning the code"""

df=min_dist_all_pat(patients)

print(df)

ids_list=pa.ids_list

patients=list(pd.read_csv(ids_list)["PatientId"])


df.to_csv(min_dist_df+"pat_test")





