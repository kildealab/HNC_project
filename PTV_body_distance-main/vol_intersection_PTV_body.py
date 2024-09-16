"""Create mask to obtain volumne intersection of PTV and body contour: 
volume = sum(slice thickness * area slice)"""
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

#import Shapely


"""Define paths"""

main_path_CT=pa.main_path_CT

main_path=pa.main_path_CBCT

plots_folder=pa.plots_folder

vol_intersection_df=pa.vol_intersection_df


def draw_all_contours(x,y,title):
    #(x,y) have to be same size

    # print("plotting ",title)
    for i in range(len(x)):
        plt.plot(x[i],y[i],"o-",label=str(i))
        plt.legend(prop={'size': 1})

        plt.gca().invert_yaxis()

        plt.savefig(plots_folder+title+str(i)+".pdf")

        plt.clf()
        plt.cla()
        plt.close()

def draw_contours(x,y,i,title):
    #(x,y) have to be same size

    # print("plotting ",title)
    # for i in range(len(x)):
    plt.plot(x[i],y[i],"o-",label=str(i))
    plt.legend(prop={'size': 1})

    plt.gca().invert_yaxis()

    plt.savefig(plots_folder+title+".pdf")


        

    plt.clf()
    plt.cla()
    plt.close()

def draw_contours_(x,y,title):
    #(x,y) have to be same size

    # print("plotting ",title)
    # for i in range(len(x)):
    plt.plot(x,y,"o-")
    plt.legend(prop={'size': 1})

    plt.gca().invert_yaxis()

    plt.savefig(plots_folder+title+".pdf")


        

    plt.clf()
    plt.cla()
    plt.close()


"""Defining functions"""

def get_CBCT_files(main_path):

	#Find the RS file

    CT_files = os.listdir(main_path)

    CT_files.sort()

    CT_files_only=[]

    for el in CT_files:

        if el[0:2]=='CT':

            CT_files_only.append(el)    

    return CT_files_only

def get_RS_file(main_path):

	#Find the RS file

    Plan_files = os.listdir(main_path)

    Plan_files.sort()

    RS_file=""

    for el in Plan_files:

        if el[0:2]=='RS':

            RS_file=el    

    return RS_file

def sort_contours(x,y,z):

	z_s,x_s,y_s=list(zip(*sorted(zip(z,x,y))))

	return x_s,y_s,z_s    


def get_index(x,y,image_position,pixel_spacing):#Based on Chloe's code
	"""For a sigle slice"""	

	x_val_pixel_index = [abs(int(round((image_position[0] - xi)/pixel_spacing[0], 0))) for xi in x]
	y_val_pixel_index = [abs(int(round((image_position[1] - yi)/pixel_spacing[1], 0))) for yi in y]

	return x_val_pixel_index,y_val_pixel_index


def get_mask(x,y,value,image_position,pixel_spacing,image_shape):#Taken from Chloe's code
    x_ind,y_ind=get_index(x,y,image_position,pixel_spacing)

    #print("x_ind",x_ind)
    #print("y_ind",y_ind)

    #draw_contours_(x_ind,y_ind,"CT_contour_new_mask")


    mask = np.zeros(image_shape, dtype=np.uint8)
	
    c = np.array(x_ind)     #Row coordinates of vertices of polygon
    r = np.array(y_ind)     #Column coordinates of vertices of polygon
    rr, cc = polygon(r, c)



    mask[rr, cc] = value       #fill area of matrix found within volume with 1s
	

    return mask

def plot_mask(mask,title):

    plt.imshow(mask)

	#plt.imshow(mask, cmap='gray', vmin=0, vmax=255)


    plt.title(title)

    plt.colorbar()

	#plt.gca().invert_yaxis()
    plt.savefig(plots_folder+"mask/"+"mask"+title+".pdf")

    plt.clf()
    plt.cla()
    plt.close()



"""Intersection"""

def intersection_(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def all_contours_single_slice(x0,y0,z0,inter):
    n=len(z0)
    m=len(inter)

    all_contours_x=[]
    all_contours_y=[]


    for j in range(m):
        cont_x=[]
        cont_y=[]
        for i in range(n):
            if z0[i]==inter[j]:
                
                cont_x=cont_x+x0[i]
                cont_y=cont_y+y0[i]
        
        all_contours_x.append(cont_x)
        all_contours_y.append(cont_y)
        
    return all_contours_x,all_contours_y

# def all_contours_inter(x0,y0,z0,inter):
#     n=len(z0)
#     m=len(inter)

#     all_contours_x=[]
#     all_contours_y=[]


#     for j in range(m):
#         cont_x=[]
#         cont_y=[]
#         for i in range(n):
#             if z0[i]==inter[j]:
                
#                 cont_x=cont_x+x0[i]
#                 cont_y=cont_y+y0[i]
        
#         all_contours_x.append(cont_x)
#         all_contours_y.append(cont_y)
        
#     return all_contours_x,all_contours_y

# def all_conts_two_structures(x0,y0,z0,x1,y1,z1):
#     inter=list(np.unique(intersection_(z0,z1)))
#     x0,y0=all_contours_single_slice(x0,y0,z0,inter)
#     x1,y1=all_contours_single_slice(x1,y1,z1,inter)
#     return x0,y0,x1,y1,inter


def all_contours_inter(x,y,z,inter):
    """Get all contours from the intersection"""
    inter_ind=[]
    for i in range(len(inter)):
        #print(inter[i])
        inter_ind=inter_ind + np.where(np.array(z)==inter[i])[0].tolist()
    
    inter_ind=sorted(inter_ind)
    x_new=[x[i] for i in inter_ind]
    y_new=[y[i] for i in inter_ind]
    z_new=[z[i] for i in inter_ind]
    
        
        
    return x_new,y_new,z_new 

def mask_per_slice(masks,z,inter,image_shape): 
    #sum all the masks per slice
    
    mask_all=[]
    for j in range(len(inter)):
        mask=np.zeros(image_shape)
        for i in range(len(z)):
            if z[i]==inter[j]:
                mask=mask+masks[i]
        mask_all.append(mask)    
            
                
    return mask_all

def anterior_half(mask_body):
    n=len(mask_body)

    y_ind,x_ind=np.where(mask_body != 0)

    if len(x_ind)>0:

        delta=(max(y_ind)-min(y_ind))/2


        anterior_half_val=int(min(y_ind)+delta)

        mask_ant_half=mask_body[0:anterior_half_val,:]

    else:
        mask_ant_half=mask_body[0:int(n/2),:]

    return mask_ant_half 


def area_out_body(mask,value_PTV,pixel_spacing):

    x_ind_PTV_out,y_ind_PTV_out=np.where(mask==value_PTV)

    area_i=np.round(pixel_spacing[0]*pixel_spacing[1]*len(x_ind_PTV_out),2)

    return area_i


def homogenize_mask(mask,value):
    
    x_ind,y_ind=np.where((mask!=0))
    for i,j in zip(x_ind,y_ind):
            mask[i,j]=value

    

    return mask


def get_mask_area(k_body,body_names,RS_path,main_path_CT,value_PTV,value_CT):

    """Read CT"""

    CT_files=get_CBCT_files(main_path_CT)

    ds=pydicom.read_file(main_path_CT+"/"+CT_files[0])

    image_position = [float(x) for x in ds.ImagePositionPatient]

    pixel_spacing = [float(x) for x in ds.PixelSpacing]

    image_shape = (ds.Rows, ds.Columns)



    print("image_position",image_position)
    print("pixel_spacing",pixel_spacing)
    print("image_shape",image_shape)

    "get contours"

    x_CT,y_CT,z_CT=fc.get_contour_all_slice(body_names[k_body],RS_path)

    x_CT,y_CT,z_CT=sort_contours(x_CT,y_CT,z_CT)

    x_PTV,y_PTV,z_PTV=fc.get_contour_all_slice("PTV_ALL",RS_path)

    x_PTV,y_PTV,z_PTV=sort_contours(x_PTV,y_PTV,z_PTV)
    print("draw contours")

    #draw_all_contours(x_CT,y_CT,"ct_cnts")



    #first component to the end of the array



    #x0,y0,x1,y1,inter1=all_conts_two_structures(x_CT,y_CT,z_CT,x_PTV,y_PTV,z_PTV)

    #print("z",inter1[slice_i])

    #print("contours x,y CT",x0[slice_i],y0[slice_i])

    #print("contours x,y PTV",x1[slice_i],y1[slice_i])

    #draw_contours_(x0[slice_i],y0[slice_i],"CBCT_contour_weird"+str(slice_i))

    #draw_contours_(x1[slice_i],y1[slice_i],"PTV_contour_weird"+str(slice_i))

    

    inter=list(np.unique(intersection_(z_CT,z_PTV)))

    print("calculating contours")

    x0,y0,z0=all_contours_inter(x_CT,y_CT,z_CT,inter)

    x1,y1,z1=all_contours_inter(x_PTV,y_PTV,z_PTV,inter)

    #print(z0)
    #print(z1)


    """Continue from here: Get all the masks and sum the ones from the same slice"""

    """
    Get CT mask -->mask_0. 
    Get PTV mask -->mask_1. 

    """

    print("calculating masks")
    

    mask_1=[get_mask(x1[slice_i],y1[slice_i],value_PTV,image_position,pixel_spacing,image_shape) for slice_i in range(len(z1))]

    mask_0=[get_mask(x0[slice_i],y0[slice_i],value_CT,image_position,pixel_spacing,image_shape) for slice_i in range(len(z0))]

    #print(z0)
    #print(inter)
    
    mask_sum_1=mask_per_slice(mask_1,z1,inter,image_shape)

    mask_sum_0=mask_per_slice(mask_0,z0,inter,image_shape)

    print("homogenizing masks")

    mask_sum_1=[homogenize_mask(mask_sum_1[i],value_PTV) for i in range(len(mask_sum_1))] 

    mask_sum_0=[homogenize_mask(mask_sum_0[i],value_CT) for i in range(len(mask_sum_0))]

    area_PTV=[area_out_body(mask_sum_i,value_PTV,pixel_spacing) for mask_sum_i in mask_sum_1]


    vol_PTV=sum(area_PTV)*0.01*3#slice thickness CT=3 mm

    print("\n number of slices")
    print(len(inter))

    print(len(area_PTV))
    # mask_1=get_mask(x1[slice_i],y1[slice_i],20,image_position,pixel_spacing,image_shape)

    # mask_0=get_mask(x0[slice_i],y0[slice_i],10,image_position,pixel_spacing,image_shape)

    print("mask sum")
    mask_sum=[mask_sum_0[i]+mask_sum_1[i] for i in range(len(inter))]
    



    """Later: take the ant half using as a reference the first CT"""

    #mask_sum_ah=[anterior_half(mask_sum_i) for mask_sum_i in mask_sum]

    """Calculate area"""

    #area=[area_out_body(mask_sum_ah_i,value_PTV,pixel_spacing) for mask_sum_ah_i in mask_sum_ah]

    area=[area_out_body(mask_sum_i,value_PTV,pixel_spacing) for mask_sum_i in mask_sum]

    print("\n\nvol_PTV",vol_PTV)

    return area,inter,vol_PTV#,mask_sum_ah
    #print("area slice "+str(slice_i),area_i)

    #plot_mask(mask_ant_half,"PTV_out_20fx area ="+str(area_i)+"mm2")


#area,inter,mask_sum=get_mask_area(5,body_names,RS_path,main_path_CT,value_PTV,value_CT)
# for i in range(len(body_names)):

#     area,inter=get_mask_area(i,body_names,RS_path,main_path_CT,value_PTV,value_CT)

#     print(body_names[i],len(area))


# print("plotting masks")

# for i in range(len(mask_sum)):
#     print(i)
#     plot_mask(mask_sum[i],"maskct+ptv, slice = "+str(inter[i])+"area ="+str(area[i])+"mm")

 
def get_areas_all_slices(body_names,RS_path,main_path_CT,value_PTV,value_CT):

    df = pd.DataFrame(columns = ['z_inter','area'])

    areas=[]
    inters=[]
    vols=[]

    for i in range(len(body_names)):

        print("\n"+body_names[i])

        area,inter,vol_PTV=get_mask_area(i,body_names,RS_path,main_path_CT,value_PTV,value_CT)
        areas.append(area)
        inters.append(inter)
        vols.append(vol_PTV)

        #if i==0:#CT


        #print("\narea",area)


    
    inters_df = pd.DataFrame(inters+areas)

    vols_df=pd.DataFrame(vols)


    #     df_CT = pd.DataFrame(np.array([inter,area]).T)

    #     df_CT.columns=['z_inter','area_'+body_names[i]]

    #     df_CT['z_inter']=df_CT['z_inter'].astype("float")

    #     if len(df)==0:
    #         df=df_CT


    #         #df.Weight.astype('int64')
    #     else:



    #         df=pd.merge(df,df_CT,how='inner', on='z_inter')

    #     print(df)

    return inters_df,vols_df



"""Run the code"""

"""Read RS""" 

RS_name="RS.BodyCT1_CBCT_1_6_11_16_20_23.dcm"

RS_path=main_path+RS_name  

#rs = pydicom.read_file(RS_path)

"""Get contours"""


#x_PTV,y_PTV,z_PTV=fc.get_contour_all_slice("PTV_ALL",RS_path)

"""Sort contours"""

#x_PTV,y_PTV,z_PTV=sort_contours(x_PTV,y_PTV,z_PTV)#This is not sorted nicely

body_names=["BODY","Body-1","Body-6","Body-11","Body-16","Body-20","Body-23"]

"""Defining PTV and values: you can assign any value different from zero. CT value must be different from PTV value """

value_PTV=20
value_CT=10

"""
Read CT file:
You do not need to find every CT file because the slice thickness is 3mm for all CTs
"""

CT_files=get_CBCT_files(main_path_CT)

print(CT_files)
ds=pydicom.read_file(main_path_CT+"/"+CT_files[0])
slice_thickness=ds.SliceThickness

print("slice_thickness CT",slice_thickness)

df_areas,vols_df=get_areas_all_slices(body_names,RS_path,main_path_CT,value_PTV,value_CT)

df_areas.to_csv(vol_intersection_df+"areas_all_fx_row_format_fullbody")

vols_df.to_csv(vol_intersection_df+"vols")

print(vols_df)