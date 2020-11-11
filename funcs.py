"""This is a module of functions for DICOM image processing Python"""
"""Authors: Aixa Andrade with collaboration of Hossein Naseri for the 3D Registration"""


import numpy as np
from matplotlib import pyplot as plt
import cv2
from scipy.interpolate import interp1d
import pydicom
import os
from scipy.ndimage import zoom
import paths as pa #import paths

"""Read plots folder path here:"""


plots_folder_path=pa.plots_folder #import plots folder path


"""The following functions were built by Aixa Andrade to draw the body contour of a DICOM image"""


def DICOMImage_GrayScale(DICOM_image_pixel):

    #The DICOM images cannot be read by OpenCV, we need to convert to grayscale (0-255)
    
    image_gray_scale = DICOM_image_pixel

    #Gets the max and min intensities from the DICOM image to reescale the DICOM image 
    max_int = image_gray_scale.max()

    min_int = image_gray_scale.min()

    #Creates a function f_int_gray that interpolates the maximum and minimun intensities to the grayscale
    f_int_gray = interp1d([min_int,max_int],[0,255])

    #Applies the function to the DICOM image

    image_gray_scale = f_int_gray(image_gray_scale)

    #for OpenCV

    image_gray_scale=image_gray_scale.astype(np.int)
    
    image_gray_scale= np.uint8(image_gray_scale)

    return image_gray_scale

def getContours(pixel_image): #function to get contours using OpenCV (based in Open CV documentation)

    image_gray_scale=DICOMImage_GrayScale(pixel_image)
    #threshold
    ret, thresh=cv2.threshold(image_gray_scale,27,30,0)#value, max_val, type

    contours, hierarchy= cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)#thershold, contour mode,contour approximation method

    return contours


def DrawContourImage(pixel_image): #function to draw contours using open cv (Based in openCV documentation)
    
    image_gray_scale=DICOMImage_GrayScale(pixel_image)
    #threshold
    ret, thresh=cv2.threshold(image_gray_scale,27,30,0)#value, max_val, type

    contours, hierarchy= cv2.findContours(thresh,cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)#thershold, contour mode,contour approximation method

    # print("Number of contours:",len(contours))

    contour_img=cv2.drawContours(imgray, contours,-1,(255,255,255),3); #-1 all contours, color, thickness=3

    return contour_img

def extract_contour(image): # extract points from a contour
    contour_x=[]
    contour_y=[]
    all_contours=getContours(image)
    for j in range(len(all_contours)):
        conts=all_contours[j]
        x_point=[conts[i][0][0] for i in range(len(conts))]
        y_point=[conts[i][0][1] for i in range(len(conts))]
        contour_x.append(x_point)
        contour_y.append(y_point)

    return contour_x,contour_y


def extract_contour_mm(image,spacing,image_pos):    # extract points from a contour
    contour_x=[]
    contour_y=[]

    all_contours=getContours(image)

    for j in range(len(all_contours)):

        conts=all_contours[j]
        x_point=[conts[i][0][0]*spacing[0] + image_pos[0] for i in range(len(conts))]
        y_point=[conts[i][0][1]*spacing[1] + image_pos[1] for i in range(len(conts))]
        contour_x.append(x_point)
        contour_y.append(y_point)

    return contour_x,contour_y


def format_contour(conts):  # extract points from a contour


    x_point=[conts[i][0][0] for i in range(len(conts))]
    y_point=[conts[i][0][1] for i in range(len(conts))]

    return x_point,y_point


def max_contour(x,y,img_modality):#Gives the maiximum contour
#Needs a lot of work
    
    #Assuming that len(x)=len(y)
    contour_len = []
    if len(x)>1: 
        
        for contour in x:
            print(contour)          
            contour_len.append(len(contour))
    else:
#       print("contour_len",len(contour))
        contour_len.append(len(x))

    #identify the index of the maximum length element/largest contour

    index_max_len=contour_len.index(max(contour_len))

#   In the CT the greatest contour around the neck corresponds to the couch

    if img_modality=='CT':
        del contour_len[index_max_len] 
        #identify the index of the maximum length element/lasrgest contour
        index_max_len=contour_len.index(max(contour_len))

    # print("index_max_len",index_max_len)
    x,y=x[index_max_len],y[index_max_len]
    return x,y

def max_contour_CBCT(image,spacing,image_pos):#Gives the maiximum contour FOR CBCT
#In CT the maximum contour includes the coush and needs more considerations
#Needs a lot of work
    
    x,y=extract_contour_mm(image,spacing,image_pos)
    #Assuming that len(x)=len(y)
    contour_len = []
    if len(x)>1: 
        
        for contour in x:
            #print(contour)         
            contour_len.append(len(contour))
    else:
#       print("contour_len",len(contour))
        contour_len.append(len(x))

    #identify the index of the maximum length element/largest contour

    index_max_len=contour_len.index(max(contour_len))



    # print("index_max_len",index_max_len)
    x,y=x[index_max_len],y[index_max_len]
    return x,y


def max_contour_CBCT_pixel(image):#Gives the maiximum contour FOR CBCT in pixels
#In CT the maximum contour includes the coush and needs more considerations
#Needs a lot of work
    
    x,y=extract_contour(image)
    #Assuming that len(x)=len(y)
    contour_len = []
    if len(x)>1: 
        
        for contour in x:
            #print(contour)         
            contour_len.append(len(contour))
    else:
#       print("contour_len",len(contour))
        contour_len.append(len(x))

    #identify the index of the maximum length element/largest contour

    index_max_len=contour_len.index(max(contour_len))



    # print("index_max_len",index_max_len)
    x,y=x[index_max_len],y[index_max_len]
    return x,y


def get_structure_index(ROI_structure, RS_filepath): #Based on Chloe Litrico previous work (2018)

    #Gets the index of the structure of a DICOM RS file

    RS_file = pydicom.read_file(RS_filepath)

    ROI_i=0
    
    #Get all ROIS 
    for i, Seq_i in enumerate(RS_file.StructureSetROISequence):

        if Seq_i.ROIName == ROI_structure: #Find the index of the selected structure

            ROI_i = i

    return ROI_i        


def get_contour_structure(ROI_structure, RS_filepath): #Based on Chloe Litrico previous work (2018)
    
    #Extracts the contour coordinates (mm) of a structure (ROI_structure) in the DICOM RS file

       
    RS_file = pydicom.read_file(RS_filepath)

    contour_coor = [] 

    ROI_i=get_structure_index(ROI_structure, RS_filepath) 
    
    #Use the index of the structure to identify the contours 
    for ROI_contour_seq in RS_file.ROIContourSequence[ROI_i].ContourSequence:

        contour_coor.append(ROI_contour_seq.ContourData) #contour coordinates
      
    contour_coor=np.array(contour_coor) #Nested arrays. The interior arrays contain the coordinates [x y z] in (mm) of each point   


    return contour_coor 

def get_contour_all_slice(ROI_structure,RS_filepath): 


    all_contours_all_slice=get_contour_structure(ROI_structure, RS_filepath)

    #Get all the contours for all the slices for a single structure
    x=[]
    y=[]
    z=[]
    for i in range(len(all_contours)):

        x_coor=all_contours_all_slice[i][0::3]
        y_coor=all_contours_all_slice[i][1::3]
        z_coor=all_contours_all_slice[i][2]

        x_coor=[float(j) for j in x_coor]
        y_coor=[float(j) for j in y_coor]
        z_coor=float(z_coor)

        x.append(x_coor)
        y.append(y_coor)
        z.append(z_coor)
    return x,y,z

def get_contour_one_slice(ROI_structure,RS_filepath,slice_num,slice_thickness):

    #Get the contour for a single slice

    x_all,y_all,z_all=get_contour_all_slice(ROI_structure,RS_filepath)

    #Get the slice in mm
    
    slice_mm=z_all[0]+slice_thickness*(slice_num) -slice_thickness
    
    ind=z_all.index(slice_mm)

    x=x_all[ind]
    y=y_all[ind]


    return x,y,slice_mm


def get_contour_isocenter_slice_CT(ROI_structure,RS_filepath,isocenter_coor):

    #Get the contour for a single slice

    x_all,y_all,z_all=get_contour_all_slice(ROI_structure,RS_filepath)

    #Get the slice in mm
       
    slice_mm=isocenter_coor[2]
    ind=z_all.index(slice_mm) #Giving a -1 +1 kind of tolerance. It is not well implemented


    x=x_all[ind]
    y=y_all[ind]


    print("slice in mm",slice_mm)

    return x,y,slice_mm

def get_contour_isocenter_slice_CBCT(ROI_structure,RS_filepath,slice_thickness):

    #Get the contour for a single slice

    x_all,y_all,z_all=get_contour_all_slice(ROI_structure,RS_filepath)

    #Get the slice in mm
    
    slice_mm=0.0 #The isocenter is at 0,0,0 in the CBCT

    ind=z_all.index(slice_mm)

    x=x_all[ind]
    y=y_all[ind]

    # print(z_all[ind])
    print("slice in mm",slice_mm)

    return x,y,slice_mm


def Registration_FR(CBCT_folderpath): #Return the Registration coordinates 
#Uses the registration coordinates, they are stored in a folderpath

    CBCT_files = os.listdir(CBCT_folderpath)
    CBCT_files.sort() #It does not work that well

#Read all CBCT slice numbers
    R_coor=[]

    for CBCT_file in CBCT_files:

        if CBCT_file[0:2] == 'RE':

            print(CBCT_file)

            ds = pydicom.read_file(CBCT_folderpath + '/' + CBCT_file)
            if ds.RegistrationSequence[1].MatrixRegistrationSequence[0].MatrixSequence[0].FrameOfReferenceTransformationMatrixType != 'RIGID':
                print ('WARNING:  NON RIGID transformation %s '%ds.RegistrationSequence[1].MatrixRegistrationSequence[0].MatrixSequence[0].FrameOfReferenceTransformationMatrixType)
            # R_coor = ds.RegistrationSequence[1].MatrixRegistrationSequence[0][0x0070,0x30a][0][[0x3006,0xc6]]
            R_coor = ds.RegistrationSequence[1].MatrixRegistrationSequence[0].MatrixSequence[0].FrameOfReferenceTransformationMatrix
            R_coor=[float(i) for i in R_coor]

        #print (CBCT_file)

    return R_coor

'''
Author: H.Naseri
following fuction reads CT ang export coordinates
'''

def get_cts(CT_files_path):
    CT_files = [os.path.join(CT_files_path, f) for f in os.listdir(CT_files_path)]
    slices = {}
    ds = None
    image_position = None
    for ct_file in CT_files:
        ds = pydicom.read_file(ct_file)

        #print(ct_file)
        

        # Check to see if it is an image file.
        # print ds.SOPClassUID
        if ds.SOPClassUID == '1.2.840.10008.5.1.4.1.1.2':
            #
            # Add the image to the slices dictionary based on its z coordinate position.
            #
            slices[ds.ImagePositionPatient[2]] = ds.pixel_array

            image_position = ds.ImagePositionPatient
            # The pixel spacing or planar resolution in the x and y directions.
            ct_pixel_spacing = ds.PixelSpacing

        else:
            pass

    # The ImagePositionPatient tag gives you the x,y,z coordinates of the center of
    # the first pixel. The slices are randomly accessed so we don't know which one
    # we have after looping through the CT slice so we will set the z position after
    # sorting the z coordinates below.

    
    #image_position = ds.ImagePositionPatient
    # print 'CT', image_position
    # Construct the z coordinate array from the image index.
    z = slices.keys()
    z = sorted(z)
    ct_z = np.array(z)

    image_position[2] = ct_z[0]

    

    # Verify z dimension spacing
    b = ct_z[1:] - ct_z[0:-1]
    # z_spacing = 2.5 # Typical spacing for our institution
    if abs(b.min() - b.max())< 1.E-6:
         z_spacing = b.max()
    else:
        print ('Error z spacing in not uniform')
        z_spacing = 0
        print (b)

    # print z_spacing

    # Append z spacing so you have x,y,z spacing for the array.
    ct_pixel_spacing.append(z_spacing)

    # Build the z ordered 3D CT dataset array.
    ct_array = np.array([slices[i] for i in z])

    # print('array shape:', ct_array.shape)    
    # plt.imshow(ct_array[50,:,:], origin='lower')
    # plt.imshow(ct_array[:,255,:], origin='lower')
    # plt.title(CT_files_path[22:])
    # plt.imshow(ct_array[:,:,255], origin='lower')
    # plt.show()
    # Now construct the coordinate arrays
    # print ct_pixel_spacing, image_position
    x = np.arange(ct_array.shape[2])*ct_pixel_spacing[0] + image_position[0]
    y = np.arange(ct_array.shape[1])*ct_pixel_spacing[1] + image_position[1]
    z = np.arange(ct_array.shape[0])*z_spacing + image_position[2]
    # print x
    # print image_position[0], image_position[1], image_position[2]
    # print ct_pixel_spacing[0], ct_pixel_spacing[1], ct_pixel_spacing[2]
    # print x, y
    # print (len(x), len(y))
    # # The coordinate of the first pixel in the numpy array for the ct is then  (x[0], y[0], z[0])
    # print('pixel_spacing: ', ct_pixel_spacing)
    return ct_array, x,y,z, ct_pixel_spacing

'''
This fucntion registers two 3D volumes
'''
def register_3d_volumes(vol1_array, ds_x,ds_y,ds_z, vol2_array, ct_x,ct_y,ct_z, cbct_ref_trans_matrix, ct_ref_trans_matrix):
    '''
    vol1_array: 3D array of dose [z,y,x] (or any other 3D volume) for example dose vol1_array[0,:,:] is 2D array 
        of first slice. vol1_array[0,0,:] is an array of pixel values of first row in 
        first slice of the dose.
    vol2_array: 3D array of CT [z,y,x]
    ds_x: ds_y,ds_z: position of x,y,z components of dose array in mm relative to isocenter. 
        For example ds_x[0] is the position of the first pixel of the vol1_array in
        x direction compared to isocenter. like if ds_x[0] = 55, ds_y[0] = -33 and ds_z[0] = 18
        the first pixel of vol1_array is in the position of [18,-33,55] mm from isocenter
    ct_x,ct_y,ct_z: osition of x,y,z components of ct array in mm relative to isocenter.
    '''
    registered_vol = None
    vol2_array = None
    '''
    cbct_t_x, cbct_t_y, cbct_t_z are transformations in x, y and z 
    between CBCT and correspoinding CT 

    ct_t_x, ct_t_y and ct_t_x are transformations in x, y and z between CT and FP1.
    for FP1 these values are zero (meaning that CT is FP1) for M1P1.

    '''
    print(cbct_ref_trans_matrix)

    print(ct_ref_trans_matrix)

    cbct_t_x = cbct_ref_trans_matrix[3]
    cbct_t_y = cbct_ref_trans_matrix[7]
    cbct_t_z = cbct_ref_trans_matrix[11]
    ct_t_x = ct_ref_trans_matrix[3]
    ct_t_y = ct_ref_trans_matrix[7]
    ct_t_z = ct_ref_trans_matrix[11]

    print("worked fine")
    '''
    ds_x, ds_y, ds_z are coordinate sytems of the CBCT images.
    Following is fixing the shift for registration between CBCT and correspoinding ct
    '''
    # commenting out the following will ignore reqiatration between CBCT and correspoinding CT
    ds_x = ds_x+cbct_t_x
    ds_y = ds_y+cbct_t_y
    ds_z = ds_z+cbct_t_z
    # testing ct1-ct2 reg
    # ds_x = ds_x+cbct_t_x+ct_t_x
    # ds_y = ds_y+cbct_t_y+ct_t_y
    # ds_z = ds_z+cbct_t_z+ct_t_z
    '''
    ct_x, ct_y, ct_z are coordinate sytems of the CT images.
    Following is fixing the shift for registration between current CT and FP1. 
    '''
    # commenting out the following will ignore reqiatration between CT and FP1
    # ct_x = ct_x+ct_t_x
    # ct_y = ct_y+ct_t_y
    # ct_z = ct_z+ct_t_z
    ct_x = ct_x
    ct_y = ct_y
    ct_z = ct_z
    x_scale = (ct_x[1]-ct_x[0])/(ds_x[1]-ds_x[0])
    y_scale = (ct_y[1]-ct_y[0])/(ds_y[1]-ds_y[0])
    z_scale = (ct_z[1]-ct_z[0])/(ds_z[1]-ds_z[0]) 
    # print (x_scale, y_scale, z_scale)
    # print ('zoom_shape', vol2_array.shape, new_v2_array.shape)
    ct_x = zoom(ct_x, x_scale)
    ct_y = zoom(ct_y, y_scale)
    ct_z = zoom(ct_z, z_scale)
    # vol2_array = zoom(vol2_array, (z_scale, y_scale, x_scale))
    vol2_array = np.zeros((len(ct_z), len(ct_y), len(ct_x)))
    # print ('zoom_shape ', vol2_array.shape, ct_z.shape, ct_y.shape, ct_x.shape)
    # print ('ds_x', np.min(ds_x), np.max(ds_x),'ds_y', np.min(ds_y), np.max(ds_y), 'ds_z', np.min(ds_z), np.max(ds_z))
    # print ('ds_x', ds_x[0], ds_x[-1],'ds_y', ds_y[0], ds_y[-1], 'ds_z', ds_z[0], ds_z[-1])
    # print ('ct_x', np.min(ct_x), np.max(ct_x), 'ct_y', np.min(ct_y), np.max(ct_y),'ct_z', np.min(ct_z), np.max(ct_z))
    # print ('ct_x', ct_x[0], ct_x[-1],'ct_y', ct_y[0], ct_y[-1], 'ct_z', ct_z[0], ct_z[-1])
    ds_0 = (np.argmin(abs(ds_x)),np.argmin(abs(ds_y)),np.argmin(abs(ds_z)))
    # find x and y cut-off for start of image
    if ds_x[0] > ct_x[0]:
        min_ct_x = np.argmin(abs(ct_x-ds_x[0]))
        min_ds_x = 0
    else:
        min_ds_x = np.argmin(abs(ds_x-ct_x[0]))
        min_ct_x = 0 
    
    if ds_x[-1] < ct_x[-1]:
        max_ct_x = np.argmin(abs(ct_x-ds_x[-1]))
        max_ds_x = len(ds_x)
    else:
        max_ct_x = np.argmin(abs(ds_x-ct_x[-1]))
        max_ds_x = len(ct_x)
    #trim max_ct
    if (max_ct_x-min_ct_x)-(max_ds_x-min_ds_x) == 1:
        max_ct_x -= 1
    elif (max_ct_x-min_ct_x)-(max_ds_x-min_ds_x) == -1:
        max_ct_x += 1
    elif (max_ct_x-min_ct_x)-(max_ds_x-min_ds_x) != 0:
        print ('Max Value Error')
    # y component of max and min cut_off 
    if ds_y[0] > ct_y[0]:
        min_ct_y = np.argmin(abs(ct_y-ds_y[0]))
        min_ds_y = 0
    else:
        min_ds_y = np.argmin(abs(ds_y-ct_y[0]))
        min_ct_y = 0
         
    if ds_y[-1] < ct_y[-1]:
        max_ct_y = np.argmin(abs(ct_y-ds_y[-1]))
        max_ds_y = len(ds_y)
    else:
        max_ct_y = np.argmin(abs(ds_y-ct_y[-1]))
        max_ds_y = len(ct_y)
    #trim max_ct
    if (max_ct_y-min_ct_y)-(max_ds_y-min_ds_y) == 1:
        max_ct_y -= 1
    elif (max_ct_y-min_ct_y)-(max_ds_y-min_ds_y) == -1:
        max_ct_y += 1
    elif (max_ct_y-min_ct_y)-(max_ds_y-min_ds_y) != 0:
        print ('Error')
    # print (ds_x, ct_x)
    # print ('-----------------------------------------')
    # print (ds_y, ct_y)
    # print (ds_z, ct_z)
    # z component of cut off
    if ds_z[0] > ct_z[0]:
        min_ct_z = np.argmin(abs(ct_z-ds_z[0]))
        min_ds_z = 0
    else:
        min_ds_z = np.argmin(abs(ds_z-ct_z[0]))
        min_ct_z = 0 
    if ds_z[-1] < ct_z[-1]:
        max_ct_z = np.argmin(abs(ct_z-ds_z[-1]))
        max_ds_z = len(ds_z)
    else:
        max_ds_z = np.argmin(abs(ds_z-ct_z[-1]))
        max_ct_z = len(ct_z)
    # print ('ds_x', min_ds_x, max_ds_x,'ds_y',  min_ds_y, max_ds_y, 'ds_z', min_ds_z, max_ds_z)
    # print ('ct_x', min_ct_x, max_ct_x, 'ct_y', min_ct_y, max_ct_y,'ct_z', min_ct_z, max_ct_z)
    # print min_ds_z, max_ds_z, min_ct_z, max_ct_z

    # print(ds_z[0], ds_z[-1], ct_z[0], ct_z[-1])
    # print len(ds_z), len(ct_z)
    registered_vol = np.zeros(vol2_array.shape) 
    # registered_vol[:] = np.nan 
    for i in range(min_ds_z,max_ds_z):   
        ds_slice = vol1_array[min_ds_z+i,min_ds_y:max_ds_y, min_ds_x:max_ds_x]
        registered_vol[min_ct_z+i,min_ct_y:min_ct_y+len(ds_slice[0]),min_ct_x:min_ct_x+len(ds_slice[1])] = ds_slice
        ds_slice = None
    return registered_vol, vol2_array, (ct_x, ct_y, ct_z)

def pad_with(vector, pad_width, iaxis, kwargs):
    pad_value = kwargs.get('padder', np.nan)
    vector[:pad_width[0]] = pad_value
    vector[-pad_width[1]:] = pad_value

def plot_it(reg_metrics,v1_reg, v2_reg, roi_z, map_name, a_ratio):
    my_cmap = plt.cm.get_cmap('jet') 
    my_cmap.set_bad(alpha=0)
    plt.imshow(v2_reg[:,int(len(reg_metrics[1])/2),:], cmap='gray', origin='lower', aspect=a_ratio)
    plt.imshow(v1_reg[:,int(len(reg_metrics[1])/2),:], cmap=my_cmap, alpha=0.5, origin='lower', aspect=a_ratio)
    plt.title(map_name)
    roi_slice = np.argmin(abs(reg_metrics[2]-roi_z))
    # print ('roi_slice', roi_slice)
    plt.plot(range(len(reg_metrics[0])), [roi_slice]*len(reg_metrics[0]), 'r-')

    plots_folder=plots_folder_path #Save the files in the plost folder
    #plt.savefig(map_name)
    plt.savefig(plots_folder+map_name)
    plt.cla()
    plt.clf()
    plt.close()

    map_name = 'slice_'+map_name
    plt.imshow(v2_reg[roi_slice,:,:], cmap='gray', origin='lower')
    plt.title(map_name)

    #plt.savefig(map_name)
    plt.savefig(plots_folder+map_name)
    plt.cla()
    plt.clf()
    plt.close()
    # print(map_name)

def register_CBCT_to_CT(CBCT_files_path, CT_files_path, roi_z, is_replan):
    vol1_array, v1_x,v1_y,v1_z, v1_pixel_spacing = get_cts(CBCT_files_path)
    vol2_array, v2_x,v2_y,v2_z, v2_pixel_spacing = get_cts(CT_files_path)
    cbct_ref_trans_matrix = Registration_FR(CBCT_files_path)

    print(cbct_ref_trans_matrix)
    map_name = 'images/%s_%s.png'%(CBCT_files_path[38:], CT_files_path[47:])
    map_name = map_name.replace('/','')
    print(map_name)
    a_ratio = v1_pixel_spacing[2]/v1_pixel_spacing[1]
    # if is_replan:
    #     print("CT_files_path",CT_files_path)
    #     ct_ref_trans_matrix = Registration_FR(CT_files_path)
    # else:
    #     ct_ref_trans_matrix = [0]*len(cbct_ref_trans_matrix) 

    #If we dont register ct1 to ct2: 

    ct_ref_trans_matrix = [0]*len(cbct_ref_trans_matrix)  

        # print('ct_ref_trans_matrix >>>', ct_ref_trans_matrix)
    v1_reg, v2_reg, reg_metrics = register_3d_volumes(vol1_array, v1_x,v1_y,v1_z, 
            vol2_array, v2_x,v2_y,v2_z, cbct_ref_trans_matrix, ct_ref_trans_matrix)
    plot_it(reg_metrics,v1_reg, v2_reg, roi_z, map_name, a_ratio)
    return v1_reg, v2_reg, reg_metrics

# def get_roi_slice(v1_reg, v2_reg, reg_metrics, roi_z):
#   roi_slice = np.argmin(abs(reg_metrics[2]-roi_z))
#   plt.plot(range(len(reg_metrics[0])), [roi_slice]*len(reg_metrics[0]), 'r-')
#   plt.savefig(map_name)
#   plt.cls()
#   plt.clf()
#   plt.close()


"""The following functions were created by Aixa Andrade"""

def obtain_slice_num_isocenter(CBCT_folderpath,isocenter):

    #It is useful to obtain the slice in which is located the real isocenter=! CBCT origin

    CBCT_files = os.listdir(CBCT_folderpath)
    CBCT_files.sort()
    z_mm=[]

    z_mm_r=[]


#Read all CBCT slice numbers


    for CBCT_file in CBCT_files:

        if CBCT_file[0:2] == 'CT':

            ds = pydicom.read_file(CBCT_folderpath + '/' + CBCT_file)
    # for slice_num in range(1,93): #field of view CBCT 93 slices

            
            #print(ds.ImagePositionPatient)

            z_mm.append(float(ds.ImagePositionPatient[2]))

            z_mm_r.append(round(float(ds.ImagePositionPatient[2]))) #because one slice did not have the exact zero



#   slice_zero=0.0 #The origin is at 0,0,0 in the CBCT

    z_mm.sort()

    z_mm_r.sort()


#   ind=z_mm_r.index(slice_zero) #index of the origin slice


    try:
        slice_zero=0.0 # Normally the origin is at 0,0,0 in the CBCT

        ind=z_mm_r.index(slice_zero)

    except ValueError:

        print("isocenter is not exactly at 0.0, it was used the slice that corresponds to 1")

        slice_zero=1 # if all the digits were rounded, then for #>0.5==1. It is used the positive part because the first slice is acquired from the feet to head, starting with negative numbers towards positive numbers. We want it to be as far as possible from the shoulders

        ind=z_mm_r.index(slice_zero) #index of the origin slice 





    slice_num_origin=ind+1 #The index goes from 0-92, slices go from 1-93

    R_coor=Registration_FR(CBCT_folderpath) #obtain the registration coordinates

    iso_slice=R_coor[11] #Z coordinate that defines the translation from the registration matrix

    delta=iso_slice-isocenter[2] #difference between origin and isocenter in CBCT

    iso_slice_shift=round(delta/ds.SliceThickness)

    slice_num_iso=slice_num_origin+iso_slice_shift #The isocenter is not at the origin. It is shifted. We found how many slices is shifted from the origin slice (0,0,0)

    # print(slice_num_iso)

    return z_mm,slice_num_iso,ds.SliceThickness



def get_shift_mid_sub(ROI_Name,RS_filepath, isocenter, slice_thickness):

    x_sub_L,y_sub_L,z_sub_L=get_contour_all_slice(ROI_Name,RS_filepath)

    coor_mid_sub=z_sub_L[-1]+ (z_sub_L[0]- z_sub_L[-1])/2 #obtaining the corrdinate of the middle subgland

    s_coor_CBCT=coor_mid_sub-isocenter[2] #CBCT shift in mm ( from isocenter to subglands)

    shift_sub_coor_CBCT=round(s_coor_CBCT/slice_thickness) # shift in slices


    return shift_sub_coor_CBCT 



def find(string,char_selec): #Based on stackoverflow answers
#Returns the index of char_select in string

    index_char=[ind for ind, char_i in enumerate(string) if char_i == char_selec]

    return index_char 

'''
Hossein:
This is the function to get the z position of the ROI slice 
'''
def get_slice_z():
    pass
def get_slice_path(CBCT_folderpath,ROI_Name,RS_filepath, isocenter):

    z_mm,slice_iso,slice_thickness=obtain_slice_num_isocenter(CBCT_folderpath,isocenter)

    shift_sub_coor_CBCT=get_shift_mid_sub(ROI_Name,RS_filepath, isocenter, slice_thickness)

    slice_final=slice_iso+shift_sub_coor_CBCT


    CBCT_files = os.listdir(CBCT_folderpath)

    name_file=CBCT_files[0] #this is just to find the names format (only valid for Aria 11)
    #the the name files selected is not CT, then take the next one. This was implemented because I got an error from the registration file
    if name_file[0:2]=='CT':

        #end_string: Get the position up to the Second point, to know the name of the file

        end_string=find(name_file,'.')[1] 
    else:

        name_file=CBCT_files[1]
    
        end_string=find(name_file,'.')[1]   

    generic_name_file=name_file[0:end_string]

    # print("position of the point:",end_string)

    # print("file name:",generic_name_file)


    selected_slice_path=str(CBCT_folderpath)+'/'+str(generic_name_file)+'.'+ str(slice_final)+'.dcm'

    # print("selected_slice_path",selected_slice_path)

    return selected_slice_path




def max_contour_CBCT_openCV(image):#Gives the maiximum contour FOR CBCT

    #Uses the format of openCV contours
    contours = getContours(image)
    #Assuming that len(x)=len(y)
    contour_len = []
    if len(contours)>1: 
        
        for cont in contours:
            #print(contour)         
            contour_len.append(len(cont))
    else:
#       print("contour_len",len(contour))
        contour_len.append(len(contours))

    #identify the index of the maximum length element/largest contour

    index_max_len=contour_len.index(max(contour_len))


    max_cont=contours[index_max_len]

    return max_cont #, contour_len



def find_CBCT_files(CBCT_main_path):

    #Obtain only the files names that corresponds to the CBCT by week

    CBCT_files = os.listdir(CBCT_main_path)

    CBCT_files.sort()

    CBCT_list=[]

    for el in CBCT_files:

        if el[0:6]=='CBCT_w':

            CBCT_list.append(el)

    return(CBCT_list)


def find_CBCT_files_general(CBCT_main_path):

    #Obtain only the files names that corresponds to the CBCT

    CBCT_files = os.listdir(CBCT_main_path)

    CBCT_files.sort()

    CBCT_list=[]

    for el in CBCT_files:

        if el[0:4]=='CBCT':

            CBCT_list.append(el)

    return(CBCT_list)

def find_RP_RS_files(main_path,PlanSetupId):

    #Obtain the names of "RP" and "RS" files given PlanSetUpId "FP1","M1P1"
    #Input variables: 
    #"PlanSetupId" --> "FP1", "M1P1"

    Plan_files = os.listdir(main_path+"/"+PlanSetupId)

    Plan_files.sort()

    RP_file=""
    RS_file=""

    for el in Plan_files:

        if el[0:2]=='RP':

            RP_file=el

        elif el[0:2]=='RS':

            RS_file=el    

    return RP_file,RS_file



    


def find_M1P1(CBCT_main_path):

    #Obtain only the file names that correspond to the CBCT files M1P1

    CBCT_files = os.listdir(CBCT_main_path)

    CBCT_files.sort()

    CBCT_list=[]

    for el in CBCT_files:

        if el=='M1P1':

            CBCT_list.append(el)

    return(CBCT_list)







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
