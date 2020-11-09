"""This script was created by Aixa Andrade to obtain the volume loss of a Region of interest of the Head and Neck  """

import os
import numpy as np
import pydicom
from matplotlib import pyplot as plt
import cv2
import funcs as fc
import pandas as pd
import paths as pa #import paths
#################################################################################################################DEfine main paths
#################################################################################################################DEfine main paths
"""Read paths here"""

images_path=pa.images_path #import image's folder path

def def_paths(patient_id,path):#patient_id on string format

	main_path=path+patient_id

	M1P1=fc.find_M1P1(main_path)

	RT_path_FP1=main_path+'/FP1/RP.'+patient_id+'.FP1 ENT.dcm'

	RT_path_M1P1=main_path+'/M1P1/RP.'+patient_id+'.M1P1 ENT.dcm'

	RS_path_FP1=main_path+'/FP1/RS.'+patient_id+'.dcm'

	RS_path_M1P1=main_path+'/M1P1/RS.'+patient_id+'.dcm'

	if len(M1P1)==0:

		print("No replan")

		paths=[RT_path_FP1,RS_path_FP1]		

	else:
		
		paths=[RT_path_FP1,RT_path_M1P1,RS_path_FP1,RS_path_M1P1]

	return paths	

			

############################################################################################# Isocenter coordinates in CT


def get_ROIName_subglands(RS_path):

	ds_RS=pydicom.read_file(RS_path)

	sub_glands_ROIName=[]
	for index, item in enumerate(ds_RS.StructureSetROISequence):

		if item.ROIName.find('Sub')!=-1:
			sub_glands_ROIName.append(item.ROIName)

	return sub_glands_ROIName	



def find_isocenter(RT_path):

	ds_RT=pydicom.read_file(RT_path) #find the isocenter in the RT file

	isocenter=ds_RT.BeamSequence[0][0x300a,0x111][0][0x300a,0x12c].value #Read the isocenter using the isocenter tag
	
	isocenter=[float(isocenter[i]) for i in range(len(isocenter))] #FRom string to float

	return isocenter		


################################################################################################## FOR LOOP


def get_vols(RT_path_FP1, RT_path_M1P1,RS_path_FP1,RS_path_M1P1,main_path):

	vols_FP1=[] #save the volumes of n_slices for each CBCT

	vols_M1P1=[]

	real_fx_number_FP1=[]

	real_fx_number_M1P1=[]

	#To find the isocenter on FP1 and M1P1 files
	isocenter_FP1=find_isocenter(RT_path_FP1)

	print("isocenter_FP1",isocenter_FP1)

	isocenter_M1P1=find_isocenter(RT_path_M1P1)


	CBCT_files=fc.find_CBCT_files_general(main_path) #Find the list of CBCT files to be used

	RS_path=''


	for i in range(len(CBCT_files)):

		print("i",i)


		############################################################
		#The following steps are useful to identify information that comes from the name of the CBCT folder 
		CBCT_files_underscore_index=fc.find(CBCT_files[i],"_") #the CBCT files have underscore that indicate the fraction number and the plan designated

		print(CBCT_files)
		print(CBCT_files_underscore_index)

		#underscore #1: PatientId
		#underscore #2: PlanReplan 10, 11
		#underscore #3: Fraction Number per Plan
		#underscore #4 (last underscore): Real Fraction Number

		CBCT_files_2nd_underscore_index=CBCT_files_underscore_index[1] 

		#check if its plan or replan

		plan_replan=CBCT_files[i][CBCT_files_2nd_underscore_index+1:CBCT_files_2nd_underscore_index+3]

		print("plan_replan:",plan_replan)


		#The next steps are to find the real fraction number

		CBCT_files_underscore_index=fc.find(CBCT_files[i],"_")

		position_last_underscore=np.max(CBCT_files_underscore_index) #the real fraction number are the two digits after the last underscore

		last_digit_CBCT_files=len(CBCT_files[i])

		real_fraction_number=CBCT_files[i][position_last_underscore+1:last_digit_CBCT_files]

		###################################################################################
		#Now read image and calculate area

		CBCT_path_i=main_path+'/'+str(CBCT_files[i]) #Individual path of each CBCT

		print('CBCT_path_i:',CBCT_path_i)

		if plan_replan=='10':

			isocenter=isocenter_FP1

			RS_path=RS_path_FP1

			sub_glands_ROIName=get_ROIName_subglands(RS_path)

		elif plan_replan=='11':	

			isocenter=isocenter_M1P1

			RS_path=RS_path_M1P1	

		sub_glands_ROIName=get_ROIName_subglands(RS_path) #get the suuuuuuubmandibular glands ROI name	

		
		#get n_slices paths
		selected_slice_path=fc.get_n_slice_path(CBCT_path_i,sub_glands_ROIName[0],RS_path, isocenter,5)
	
		print("selected slice path:",selected_slice_path)
		#Sum the individual volumes for each slice. vol=area_i(cm2)*SliceThickness (cm)

		vol=0

		for selected_slice_path_i in selected_slice_path:

			ds_CBCT_i = pydicom.read_file(selected_slice_path_i) #Read that slice

			SliceThickness=ds_CBCT_i.SliceThickness

			pixel_CBCT_i=ds_CBCT_i.pixel_array

			pixel_spacing_i=float(ds_CBCT_i.PixelSpacing[0]) #Get pixel spacing

			pixel_area_i=pixel_spacing_i*pixel_spacing_i #Get area of each pixel

			cont_i=fc.max_contour_CBCT_openCV(pixel_CBCT_i)  #Get the body contour

			area_i = cv2.contourArea(cont_i) #Get the area of the body contour

			area_i_cm2=area_i*pixel_area_i*.01 #Pass to cm2

			vol=vol +area_i_cm2*SliceThickness*.1 #multiply by slice thickness and sum the 5 slices
			####################################################################

		if plan_replan=='10':

			vols_FP1.append(vol) #fill the vol array

			real_fx_number_FP1.append(int(real_fraction_number))

		elif plan_replan=='11':

			vols_M1P1.append(vol) #fill the vol array

			real_fx_number_M1P1.append(int(real_fraction_number))


	#sort by fraction number (they are desprganized)	
	real_fx_number_FP1, vols_FP1 = list(zip(*sorted(zip(real_fx_number_FP1, vols_FP1))))	

	real_fx_number_M1P1, vols_M1P1 = list(zip(*sorted(zip(real_fx_number_M1P1, vols_M1P1))))



	real_fx_number=real_fx_number_FP1+real_fx_number_M1P1

	vols=vols_FP1+vols_M1P1


	return CBCT_files,real_fx_number,vols


def get_vols_no_replan(RT_path_FP1,RS_path_FP1,main_path):


	real_fx_number_FP1=[]

	vols_FP1=[] #save the volumes of n_slices for each CBCT

	#To find the isocenter on FP1 and M1P1 files
	isocenter_FP1=find_isocenter(RT_path_FP1)

	CBCT_files=fc.find_CBCT_files_general(main_path) #Find the list of CBCT files to be used


	for i in range(len(CBCT_files)):


		############################################################
		#The following steps are useful to identify information that comes from the name of the CBCT folder 
		CBCT_files_underscore_index=fc.find(CBCT_files[i],"_") #the CBCT files have underscore that indicate the fraction number and the plan designated

		print(CBCT_files)
		print(CBCT_files_underscore_index)

		#underscore #1: PatientId
		#underscore #2: PlanReplan 10, 11
		#underscore #3: Fraction Number per Plan
		#underscore #4 (last underscore): Real Fraction Number

		CBCT_files_2nd_underscore_index=CBCT_files_underscore_index[1] 

		#check if its plan or replan

		plan_replan=CBCT_files[i][CBCT_files_2nd_underscore_index+1:CBCT_files_2nd_underscore_index+3]

		print("plan_replan:",plan_replan)


		#The next steps are to find the real fraction number

		position_last_underscore=np.max(CBCT_files_underscore_index) #the real fraction number are the two digits after the last underscore

		last_digit_CBCT_files=len(CBCT_files[i])

		real_fraction_number=CBCT_files[i][position_last_underscore+1:last_digit_CBCT_files]

		real_fx_number_FP1.append(int(real_fraction_number))

		###################################################################################
		#Now read image and calculate area

		CBCT_path_i=main_path+'/'+str(CBCT_files[i]) #Individual path of each CBCT

		print('CBCT_path_i:',CBCT_path_i)


		sub_glands_ROIName=get_ROIName_subglands(RS_path_FP1) #get the suuuuuuubmandibular glands ROI name	

		
		#get n_slices paths
		selected_slice_path=fc.get_n_slice_path(CBCT_path_i,sub_glands_ROIName[0],RS_path_FP1, isocenter_FP1,5)
	
		print("selected slice path:",selected_slice_path)
		#Sum the individual volumes for each slice. vol=area_i(cm2)*SliceThickness (cm)

		vol=0

		for selected_slice_path_i in selected_slice_path:

			ds_CBCT_i = pydicom.read_file(selected_slice_path_i) #Read that slice

			SliceThickness=ds_CBCT_i.SliceThickness

			pixel_CBCT_i=ds_CBCT_i.pixel_array

			pixel_spacing_i=float(ds_CBCT_i.PixelSpacing[0]) #Get pixel spacing

			pixel_area_i=pixel_spacing_i*pixel_spacing_i #Get area of each pixel

			cont_i=fc.max_contour_CBCT_openCV(pixel_CBCT_i)  #Get the body contour

			area_i = cv2.contourArea(cont_i) #Get the area of the body contour

			area_i_cm2=area_i*pixel_area_i*.01 #Pass to cm2

			vol=vol +area_i_cm2*SliceThickness*.1 #multiply by slice thickness and sum the 5 slices
			####################################################################

		vols_FP1.append(vol) #fill the areas array

	print("vols_FP1",vols_FP1)
	print("real fx number FP1",real_fx_number_FP1)




	#sort by fraction number (they are desprganized)	
	real_fx_number_FP1, vols_FP1 = list(zip(*sorted(zip(real_fx_number_FP1, vols_FP1))))	



	return CBCT_files,real_fx_number_FP1, vols_FP1


def normal_plot(x,y,x_label,y_label,title):

	plt.plot(x,y,'bo-')

	plt.xlabel(x_label)

	plt.ylabel(y_label)

	plt.title(title)

	plt.savefig(title+'.pdf')

	plt.close()


def get_vols_replanned_or_not(patient_id,path):

	main_path=path+patient_id

	paths=def_paths(patient_id,path)

	if len(paths)==2:

		RT_path_FP1,RS_path_FP1=paths

		CBCT_files,real_fx_number,vols=get_vols_no_replan(RT_path_FP1,RS_path_FP1,main_path)

	else:

		RT_path_FP1,RT_path_M1P1,RS_path_FP1,RS_path_M1P1=paths
		
		CBCT_files,real_fx_number,vols=get_vols(RT_path_FP1, RT_path_M1P1,RS_path_FP1,RS_path_M1P1,main_path)


	return CBCT_files,real_fx_number,vols	



#CBCT_files,real_fx_number,vols=get_vols_replanned_or_not(patient_id,path)

#normal_plot(real_fx_number,vols,"fx number","vol (cm3)","plot volume loss patient A")	

"""
References: 


Chloe, Litrico; John, Kildea; Haley, Patrick. 2018. Radiomics for Prostate Cancer: https://www.authorea.com/users/232204/articles/296788-radiomics-for-prostate-cancer


DICOM INNOLITICS: https://dicom.innolitics.com/ciods

DICOM STANDARD: https://www.dicomstandard.org/current

OpenCV: https://opencv.org/about/

OS: https://docs.python.org/3/library/os.html

PyDICOM: https://pydicom.github.io/

SCIPY: https://www.scipy.org/docs.html


Stackoverflow: https://stackoverflow.com/questions

"""