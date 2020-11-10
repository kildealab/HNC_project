"""This is a script that process DICOM images and returns the area in a Region of Interest of the neck"""
"""Authors: """
"""Aixa Andrade: General structure of the code and Registration in z axis"""
"""Hossein Naseri: Implementation of Registration in 3D"""
import os
import numpy as np
import pydicom
from matplotlib import pyplot as plt
import cv2
import funcs as fc
import pandas as pd
import paths as pa #import paths
#################################################################################################################DEfine main paths
"""Read paths here"""

images_path=pa.images_path #import image's folder path
plots_folder=pa.plots_folder #import plot's folder path
ids_list=pa.ids_list


def def_paths(patient_id,path):#patient_id on string format

    main_path=path+patient_id

    M1P1=fc.find_M1P1(main_path)

    FP1_RP_file,FP1_RS_file=fc.find_RP_RS_files(main_path,"FP1")

    RT_path_FP1=main_path+'/FP1/'+FP1_RP_file
    
    RS_path_FP1=main_path+'/FP1/'+FP1_RS_file

    if len(M1P1)==0:

        print("NO REPLAN")

        paths=[RT_path_FP1,RS_path_FP1]     

    else:


        M1P1_RP_file,M1P1_RS_file=fc.find_RP_RS_files(main_path,"M1P1")

        RT_path_M1P1=main_path+'/M1P1/'+M1P1_RP_file

        RS_path_M1P1=main_path+'/M1P1/'+M1P1_RS_file
        
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

def find_right_index_submand(sub_glands_ROIName):

    """Sometimes there are 4 structures for the submandibular glands:
        Example: SubmndSalv_R, SubmndSalv_L, SubmndSalv_R_NOS and SubmndSalv_L_NOS
        "NOS" do not have contours. This is a very special case but it can crash the code if not implemented"""

    index_N=[]
    for i in range(len(sub_glands_ROIName)):

        str_number=fc.find(sub_glands_ROIName[i],"N")

        if not str_number:
    
            index_N.append("good_structure") #The one without NOS
        else:
            index_N.append("weird_structure") #The one with NOS"

    return index_N.index("good_structure")



     


def find_isocenter(RT_path):

    ds_RT=pydicom.read_file(RT_path) #find the isocenter in the RT file

    isocenter=ds_RT.BeamSequence[0][0x300a,0x111][0][0x300a,0x12c].value #Read the isocenter using the isocenter tag
    
    isocenter=[float(isocenter[i]) for i in range(len(isocenter))] #FRom string to float

    return isocenter        


################################################################################################## FOR LOOP


def get_area(RT_path_FP1, RT_path_M1P1,RS_path_FP1,RS_path_M1P1,main_path):

    areas_FP1=[]

    areas_M1P1=[]

    real_fx_number_FP1=[]

    real_fx_number_M1P1=[]

    #To find the isocenter on FP1 and M1P1 files
    isocenter_FP1=find_isocenter(RT_path_FP1)
    isocenter_M1P1=find_isocenter(RT_path_M1P1)
    # print("isocenters_FP1, M1P1",isocenter_FP1, isocenter_M1P1)

    CBCT_files=fc.find_CBCT_files_general(main_path) #Find the list of CBCT files to be used

    RS_path=''


    for i in range(len(CBCT_files)):

        ############################################################
        #The following steps are useful to identify information that comes from the name of the CBCT folder 
        CBCT_files_underscore_index=fc.find(CBCT_files[i],"_") 
        #the CBCT files have underscore that indicate the fraction number and the plan designated
        # print(i, CBCT_files[i])
        # print(CBCT_files_underscore_index)
        #underscore #1: PatientId
        #underscore #2: PlanReplan 10, 11
        #underscore #3: Fraction Number per Plan
        #underscore #4 (last underscore): Real Fraction Number

        CBCT_files_2nd_underscore_index=CBCT_files_underscore_index[1] 

        #check if its plan or replan

        plan_replan=CBCT_files[i][CBCT_files_2nd_underscore_index+1:CBCT_files_2nd_underscore_index+3]
        # print("plan_replan:",plan_replan)


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
            is_replan = False
        elif plan_replan=='11': 
            isocenter=isocenter_M1P1
            RS_path=RS_path_M1P1
            is_replan = True    
        CT_path = os.path.split(RS_path)[0]
        # print('>>> CT_path', CT_path)
        sub_glands_ROIName=get_ROIName_subglands(RS_path) #get the suuuuuuubmandibular glands ROI name  
        ind_sub=find_right_index_submand(sub_glands_ROIName)
        # print ('>> isocenter is :', isocenter)
        # map_name = '%s_%s.png'%(CBCT_path_i[38:], CT_path[47:])
        # map_name = map_name.replace('/','')
        all_contours_ROI = fc.get_contour_structure(sub_glands_ROIName[ind_sub], RS_path)
        z_smg = np.mean([float(roi[2]) for roi in all_contours_ROI])
        # print(map_name)

        '''
        Aixa's Method (Registration in z axis)
        '''
        # # print (all_contours_ROI)
        # # select slice of the CBCT, considering isocenter and shift from the isocenter
        # selected_slice_path_i=fc.get_slice_path(CBCT_path_i,sub_glands_ROIName[0],RS_path, isocenter) 
        # print(selected_slice_path_i)
        # ds_CBCT_i = pydicom.read_file(selected_slice_path_i) #Read that slice
        # pixel_CBCT_i=ds_CBCT_i.pixel_array
        # pixel_spacing_x=float(ds_CBCT_i.PixelSpacing[0]) #Get pixel spacing
        # pixel_spacing_y=float(ds_CBCT_i.PixelSpacing[0]) #Get pixel spacing
        
        '''
        HOSSEIN's Method (Registration in 3D)
        '''
        reg_cbct_array, reg_ct_array, reg_metrics = fc.register_CBCT_to_CT(CBCT_path_i, CT_path, z_smg, is_replan)
        roi_slice = np.argmin(abs(reg_metrics[2]-z_smg))
        # print ('z_smg', z_smg, roi_slice)
        pixel_CBCT_i = reg_cbct_array[roi_slice,:,:]
        pixel_spacing_x = reg_metrics[0][1]-reg_metrics[0][0]
        pixel_spacing_y = reg_metrics[1][1]-reg_metrics[1][0]
        '''
        END
        '''

        pixel_area_i=pixel_spacing_x*pixel_spacing_y #Get area of each pixel
        cont_i=fc.max_contour_CBCT_openCV(pixel_CBCT_i)  #Get the body contour
        area_i = cv2.contourArea(cont_i) #Get the area of the body contour
        area_i_cm2=area_i*pixel_area_i*.01 #Pass to cm2

        ####################################################################

        if plan_replan=='10':
            areas_FP1.append(area_i_cm2) #fill the areas array
            real_fx_number_FP1.append(int(real_fraction_number))
        elif plan_replan=='11':
            areas_M1P1.append(area_i_cm2) #fill the areas array
            real_fx_number_M1P1.append(int(real_fraction_number))
        ####################################################################
        #Plot (just to get an idea)

        cnt_name = 'cnt_%s_%s_%i.png'%(CBCT_path_i[38:], CT_path[47:], i)
        cnt_name = cnt_name.replace('/','')
        x_i,y_i = fc.format_contour(cont_i)
        plt.plot(x_i,y_i,'r')
        plt.imshow(pixel_CBCT_i, cmap='gray', vmin=0, vmax=255)
        plt.title(cnt_name)
    
        plt.savefig(plots_folder+cnt_name)
        plt.clf()
        plt.cla()
        plt.close()
        # plt_contour_series(i,cont_i,pixel_CBCT_i,CBCT_files) #uses the plot function defined above
# 


    #sort by fraction number (they are desprganized)    
    real_fx_number_FP1, areas_FP1 = list(zip(*sorted(zip(real_fx_number_FP1, areas_FP1))))  
    real_fx_number_M1P1, areas_M1P1 = list(zip(*sorted(zip(real_fx_number_M1P1, areas_M1P1))))

    Plan_Replan=['10']*len(real_fx_number_FP1)+['11']*len(real_fx_number_M1P1)
    real_fx_number=real_fx_number_FP1+real_fx_number_M1P1
    areas=areas_FP1+areas_M1P1

    # print ('-------------------------------')
    # # print(CBCT_files)
    # print(real_fx_number)
    # print(areas)
    return CBCT_files,real_fx_number,areas,Plan_Replan



def get_area_no_replan(RT_path_FP1,RS_path_FP1,main_path):

    areas_FP1=[]

    real_fx_number_FP1=[]

    #To find the isocenter on FP1 and M1P1 files
    isocenter_FP1=find_isocenter(RT_path_FP1)

    CBCT_files=fc.find_CBCT_files_general(main_path) #Find the list of CBCT files to be used


    for i in range(len(CBCT_files)):

        # print("i",i)


        ############################################################
        #The following steps are useful to identify information that comes from the name of the CBCT folder 
        CBCT_files_underscore_index=fc.find(CBCT_files[i],"_") #the CBCT files have underscore that indicate the fraction number and the plan designated

        # print(CBCT_files)
        # print(CBCT_files_underscore_index)

        #underscore #1: PatientId
        #underscore #2: PlanReplan 10, 11
        #underscore #3: Fraction Number per Plan
        #underscore #4 (last underscore): Real Fraction Number

        CBCT_files_2nd_underscore_index=CBCT_files_underscore_index[1] 

        #check if its plan or replan

        plan_replan=CBCT_files[i][CBCT_files_2nd_underscore_index+1:CBCT_files_2nd_underscore_index+3]

        # print("plan_replan:",plan_replan)


        #The next steps are to find the real fraction number

        position_last_underscore=np.max(CBCT_files_underscore_index) #the real fraction number are the two digits after the last underscore

        last_digit_CBCT_files=len(CBCT_files[i])

        real_fraction_number=CBCT_files[i][position_last_underscore+1:last_digit_CBCT_files]

        ###################################################################################
       ###################################################################################
        #Now read image and calculate area
        CBCT_path_i=main_path+'/'+str(CBCT_files[i]) #Individual path of each CBCT
        # print('CBCT_path_i:',CBCT_path_i)
        if plan_replan=='10':
            isocenter=isocenter_FP1
            RS_path=RS_path_FP1
            sub_glands_ROIName=get_ROIName_subglands(RS_path)
            is_replan = False
        elif plan_replan=='11': 
            isocenter=isocenter_M1P1
            RS_path=RS_path_M1P1
            is_replan = True    
        CT_path = os.path.split(RS_path)[0]
        # print('>>> CT_path', CT_path)
        sub_glands_ROIName=get_ROIName_subglands(RS_path) #get the suuuuuuubmandibular glands ROI name  
        ind_sub=find_right_index_submand(sub_glands_ROIName)
        # print ('>> isocenter is :', isocenter)
        # map_name = '%s_%s.png'%(CBCT_path_i[38:], CT_path[47:])
        # map_name = map_name.replace('/','')
        all_contours_ROI = fc.get_contour_structure(sub_glands_ROIName[ind_sub], RS_path)
        z_smg = np.mean([float(roi[2]) for roi in all_contours_ROI])
        # print(map_name)

        '''
        Aixa's Method (Registration in z axis)
        '''
        # # print (all_contours_ROI)
        # # select slice of the CBCT, considering isocenter and shift from the isocenter
        # selected_slice_path_i=fc.get_slice_path(CBCT_path_i,sub_glands_ROIName[0],RS_path, isocenter) 
        # print(selected_slice_path_i)
        # ds_CBCT_i = pydicom.read_file(selected_slice_path_i) #Read that slice
        # pixel_CBCT_i=ds_CBCT_i.pixel_array
        # pixel_spacing_x=float(ds_CBCT_i.PixelSpacing[0]) #Get pixel spacing
        # pixel_spacing_y=float(ds_CBCT_i.PixelSpacing[0]) #Get pixel spacing
        
        '''
        HOSSEIN's Method (Registration in 3D)
        '''
        reg_cbct_array, reg_ct_array, reg_metrics = fc.register_CBCT_to_CT(CBCT_path_i, CT_path, z_smg, is_replan)
        roi_slice = np.argmin(abs(reg_metrics[2]-z_smg))
        # print ('z_smg', z_smg, roi_slice)
        pixel_CBCT_i = reg_cbct_array[roi_slice,:,:]
        pixel_spacing_x = reg_metrics[0][1]-reg_metrics[0][0]
        pixel_spacing_y = reg_metrics[1][1]-reg_metrics[1][0]
        '''
        END
        '''

        pixel_area_i=pixel_spacing_x*pixel_spacing_y #Get area of each pixel
        cont_i=fc.max_contour_CBCT_openCV(pixel_CBCT_i)  #Get the body contour
        area_i = cv2.contourArea(cont_i) #Get the area of the body contour
        area_i_cm2=area_i*pixel_area_i*.01 #Pass to cm2

        areas_FP1.append(area_i_cm2) #fill the areas array

        real_fx_number_FP1.append(int(real_fraction_number))

    #sort by fraction number (they are desprganized)    
    real_fx_number_FP1, areas_FP1 = list(zip(*sorted(zip(real_fx_number_FP1, areas_FP1))))  
    Plan_Replan=['10']*len(real_fx_number_FP1)

    return CBCT_files,real_fx_number_FP1, areas_FP1,Plan_Replan


def normal_plot(x,y,x_label,y_label,title):

    
    plt.plot(x,y,'bo-')

    plt.xlabel(x_label)

    plt.ylabel(y_label)

    plt.title(title)

    

    plt.savefig(plots_folder+title+'.pdf')
    plt.clf()
    plt.cla()
    plt.close()

def get_areas_replanned_or_not(patient_id,path):

    main_path=path+patient_id

    paths=def_paths(patient_id,path)

    if len(paths)==2:
        print ('NO REPLAN')
        RT_path_FP1,RS_path_FP1=paths
        CBCT_files,real_fx_number,areas,Plan_Replan=get_area_no_replan(RT_path_FP1,RS_path_FP1,main_path)
    else:
        print ('WITH REPLAN')
        RT_path_FP1,RT_path_M1P1,RS_path_FP1,RS_path_M1P1=paths
        CBCT_files,real_fx_number,areas,Plan_Replan=get_area(RT_path_FP1, RT_path_M1P1,RS_path_FP1,RS_path_M1P1,main_path)

    return CBCT_files,real_fx_number,areas,Plan_Replan  

def get_areas_all_patients(patients_list,path):
    #function to generate a dataframe with the area obtained from a list of patients

    column_names={'PatientId','real_fx_number','area','Plan_Replan'}

    data = pd.DataFrame(columns = column_names)

    for patient_id in patients_list:


        CBCT_files,real_fx_number,areas,Plan_Replan=get_areas_replanned_or_not(patient_id,path) #The CBCT files are not in order


        df_dictionary = {'PatientId':[patient_id]*len(real_fx_number),'real_fx_number':real_fx_number,'area':areas,'Plan_Replan':Plan_Replan}

        df= pd.DataFrame(df_dictionary)

        data = pd.concat([data, df])
        
        
        #data.to_csv('dataframe_area_loss_%s'%patient_id)
        normal_plot(real_fx_number,areas,"fx number","area (cm2)","plot area loss patient_%s"%patient_id)
    
    data.to_csv('dataframe_area_loss')
    return data   


def put_zeros_left(data,col):
    #If the PatientId is in number format lose the zeros on the left             
    for i in range(len(data)):


        len_str=len(str(data.loc[i,str(col)]))
        if len_str<7:
            diff=7-len_str
            more_zeros='0'*diff
            data.loc[i,str(col)]=more_zeros+str(data.loc[i,str(col)])

        else:
            data.loc[i,str(col)]=str(data.loc[i,str(col)])


    return data

def format_patient_id_list(patients_file_name): 
    #python does not read zeros on the left, also, it does not identify the elemnts as strings
 
    patient_id_list=pd.read_csv(patients_file_name)
    #patient_id_list=put_zeros_left(patient_id_list,'PatientId')    
    patient_id_list=patient_id_list['PatientId'].values.tolist()
    patient_id_list=[str(i) for i in patient_id_list]  

    return(patient_id_list)



def main():
    path=images_path #Define the images path
   
    patients=format_patient_id_list(ids_list) #insert a csv file with a column of "PatientId"
#    patients.remove(patient_id) #If you want to remove a patient_id from the list
    get_areas_all_patients(patients,path)


if __name__=="__main__":
    main()






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