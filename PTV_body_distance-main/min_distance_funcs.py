"""This script is about including all body contours 
Implemented negative min distance for points outside polygon
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






"""Defining functions"""


def sort_contours(x,y,z):

	z_s,x_s,y_s=list(zip(*sorted(zip(z,x,y))))

	return x_s,y_s,z_s    


"""Draw contours"""

def draw_contours(x,y,z,title):
	#(x,y) have to be same size

	print("plotting ",title)
	for i in range(len(x)):
		plt.plot(x[i],y[i],label=str(z[i]))
		plt.legend(prop={'size': 1})

		plt.gca().invert_yaxis()

		plt.savefig(plots_folder+title+".pdf")


		

	plt.clf()
	plt.cla()
	plt.close()

def draw_contours_range(x1,y1,x2,y2,z1,z2,start,stop,step, title):
	#(x,y) have to be same size
	"""Example: 
	CT-->1, CBCT-->2

	CBCT-->1, PTV-->2
	"""

	print("plotting body contours")

	for i in range(start,stop,step):
	
		plt.plot(x1[i],y1[i],"--",label="RS1 z="+str(z1[i])+" mm")

		plt.plot(x2[i],y2[i],"-",label="RS2 z="+str(z2[i])+"mm")



		plt.xlabel("x (mm)")

		plt.ylabel("y (mm)")

		plt.title(title)

		plt.legend(prop={'size': 1.5},bbox_to_anchor=(1.05, 1))


	plt.gca().invert_yaxis()
	plt.savefig(plots_folder+title+str(i)+".pdf")  



	plt.clf()
	plt.cla()
	plt.close()

def plot_min_distance(x1,y1,x0,y0,min_dist_index_1,min_dist_index_0,title):
	#(x,y) have to be same size
	"""Example: 


	CBCT-->1, PTV-->0
	"""

	print("plotting body contours")


	
	plt.plot(x1,y1,"--+",label="CT ")

	plt.plot(x0,y0,"-",label="PTV (anterior half) ")

	x=[x0[min_dist_index_0],x1[min_dist_index_1]]

	y=[y0[min_dist_index_0],y1[min_dist_index_1]]



	plt.plot(x,y,"ro-",label="min distance")


	plt.xlabel("x (mm)")

	plt.ylabel("y (mm)")

	plt.title(title)

	plt.legend(prop={'size': 12})#)


	plt.gca().invert_yaxis()
	plt.savefig(plots_folder+"min_dist/"+title+".pdf")  



	plt.clf()
	plt.cla()
	plt.close()

def plot_poly(polygon):
	x,y = polygon.exterior.xy
	plt.plot(x,y)
	plt.gca().invert_yaxis()
	plt.savefig(plots_folder+"min_dist/plot_poly.pdf")  
	plt.clf()
	plt.cla()
	plt.close()


"""Intersection (Not working because only takes half of the PTV contours- need to fix it"""

def intersection(lst1, lst2): #Taken from geeksforgeeks.com

	lst1=[float(x) for x in lst1]

	lst2=[float(x) for x in lst2]

	inter=list(set(lst1).intersection(lst2))

	inter=sorted(inter)

	ind_1 = [lst1.index(x) for x in inter]
	
	ind_2 = [lst2.index(x) for x in inter]

	return  inter,ind_1,ind_2


def inter_RS_multi_cnt_slice(x1,y1,z1,x2,y2,z2):
	"""
	Assumes than R1 structure is greater than R2 and all the values z2 are in z1 
	Example
	R1-->CT
	R2-->PTV (CBCT has more than one contour per slice)
	It will return only common slices between R1,R2
	"""	

	z_inter,ind_1,ind_2=intersection(z1,z2)

	x_1_int=[x1[i] for i in ind_1]
	y_1_int=[y1[i] for i in ind_1]
	z_1_int=[z1[i] for i in ind_1]
	""""
	Even if 2 structures have common slices(z1,z2), they can have more than one contour for slice
	-->each structure own z set
	"""

	x_2_int=[x2[i] for i in ind_2]
	y_2_int=[y2[i] for i in ind_2]
	z_2_int=[z2[i] for i in ind_2]

	return x_1_int,y_1_int,z_1_int,x_2_int,y_2_int,z_2_int


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

def all_conts_two_structures(x0,y0,z0,x1,y1,z1):
    inter=list(np.unique(intersection_(z0,z1)))
    x0,y0=all_contours_single_slice(x0,y0,z0,inter)
    x1,y1=all_contours_single_slice(x1,y1,z1,inter)
    return x0,y0,x1,y1,inter


def anterior_half_cont(y):

	#print("\ncalculating anterior threshold")
	#print("y0",y)
	#single slice
	#print("min",min(y))
	delta=(max(y)-min(y))/2

	#print("delta",delta)
	anterior_half_y=int(min(y)+delta)
	return anterior_half_y

	

def get_points_ant_half(x,y,y0):
	#y0=anterior_half_cont(y)
	y_ant_half=[yi for yi in y if yi<y0]
	x_ant_half=[xi for xi,yi in zip(x,y) if yi<y0]
    
	return x_ant_half,y_ant_half

def get_ant_half_all_slices(y):

	anterior_half_list=[]

	for i in range(len(y)):

		anterior_half_list.append(anterior_half_cont(y[i]))

	return anterior_half_list


def apply_ah(x,y,anterior_half_y0):
	#single slice


	y_ah = [y_i for y_i in y if y_i < anterior_half_y0]
	x_ah=x[0:len(y_ah)]

	return x_ah,y_ah

def apply_ah_all_slices(x,y,anterior_half_list):
	#input list of lists

	x_ah=[]
	y_ah=[]

	for i in range(len(x)):

		x_ah_i,y_ah_i=apply_ah(x[i],y[i],anterior_half_list[i])

		x_ah.append(x_ah_i)
		y_ah.append(y_ah_i)

	return x_ah,y_ah





def max_cnt_CT(x,y,z):
    
    max_cnt_indexes=[]
    
    x_new=[]
    
    y_new=[]
    
    z_unique=np.unique(z)
    
    for i in range(len(z_unique)):
    
        ind=np.where(np.array(z)==z_unique[i])[0].tolist() #index of the slices with same z position
    
        cnt_len=[len(x[i]) for i in ind] #lenght of the x coordinates of the slices with the same z postion
    
        max_len_ind=ind[cnt_len.index(max(cnt_len))]#select the index of the the largest contour
    
        max_cnt_indexes.append(max_len_ind)
    
    """end for loop"""
    
    x_new=[x[i] for i in max_cnt_indexes]
    
    y_new=[y[i] for i in max_cnt_indexes]
    
    
    return x_new,y_new,z_unique
    
    
     


"""Minimum distance"""

def d(x1,y1,x0,y0):
    return np.sqrt((x1-x0)**2 + (y1-y0)**2)

def distance(x1,y1,x0,y0):

    """CT-->1, PTV-->0"""

    d_min_one_all=[]
    
    ind_y_d_min_one_all=[]
    
    #body contour polygon

    polygon=[tuple([xi,yi]) for xi,yi in zip(x1,y1)]

    line = geometry.LineString(polygon)
    
    polygon = geometry.Polygon(line)

    


    #iterate for all the PTV point
    
    for x0_j,y0_j in zip(x0,y0):

        d_one_point_vs_all_points=[]


        point = geometry.Point(np.round(x0_j,1), np.round(y0_j,1))
        
        
        if polygon.contains(point):
            
            #print("inside")
            
            d_one_point_vs_all_points=[d(x1_i,y1_i,x0_j,y0_j) for x1_i,y1_i in zip(x1,y1)]
        
            if len(d_one_point_vs_all_points)>0:

                min_one_all=min(d_one_point_vs_all_points)                 

        else:
            #print("\outside")

            d_one_point_vs_all_points=[-1*d(x1_i,y1_i,x0_j,y0_j) for x1_i,y1_i in zip(x1,y1)]
            
            if len(d_one_point_vs_all_points)>0:
            
                min_one_all=max(d_one_point_vs_all_points)
                    

        d_min_one_all.append(min_one_all)
        
        ind_y_d_min_one_all.append(d_one_point_vs_all_points.index(min_one_all))


    if len(d_min_one_all)>0:
                
        min_dist=min(d_min_one_all)
        min_dist_index_0=d_min_one_all.index(min_dist)
        min_dist_index_1=ind_y_d_min_one_all[min_dist_index_0]

    else:
        min_dist="NA"
        min_dist_index_0="NA"
        min_dist_index_1="NA"
        
    print("\n min dist",min_dist)


    return min_dist,min_dist_index_1,min_dist_index_0

def min_distance_all_slices(x1,y1,x0,y0,inter,title):
	"""Example CT-->1, PTV-->0
	Technically, after intersecting RS structures,z1=z0"""

	min_dist_arr=[]
	min_dist_index_1_arr=[]
	min_dist_index_0_arr=[]



	"Calculating min distance for all slices \n"

	#for k in range(0,2): #Problems in slices 24,25
	for k in range(len(inter)):#iteration over all slices

		print("\nslice",k)
		print("z",inter[k])


		# x0_=x0[k]
		# y0_=y0[k]

		y_ant_half_threshold=anterior_half_cont(y1[k])

		x0_,y0_=get_points_ant_half(x0[k],y0[k],y_ant_half_threshold) #anterior half points only

		#min_dist,min_dist_index_1,min_dist_index_0=distance(x1[k],y1[k],x0[k],y0[k])
		min_dist,min_dist_index_1,min_dist_index_0=distance(x1[k],y1[k],x0_,y0_)

		min_dist_arr.append(min_dist)
		min_dist_index_1_arr.append(min_dist_index_1)
		min_dist_index_0_arr.append(min_dist_index_0)

		# print("\n min_dist,min_dist_index 1,0:")
		# print(min_dist,min_dist,min_dist_index_1,min_dist_index_0)
		# plt.plot(x0[k],y0[k],"-")
		# plt.plot(x0_,y0_,"-")

		# plot_min_distance(x1[k],y1[k],x0_,y0_,min_dist_index_1,min_dist_index_0,"slice "+str(inter[k])+title+" "+str(np.round(min_dist,2))+" mm")

		# plt.savefig(plots_folder+"PTV"+str(k)+".pdf")

		# plt.clf()
		# plt.cla()
		# plt.close()

	"""start of comment"""
	print(inter,min_dist_arr)

	df = pd.DataFrame(np.array([inter,min_dist_arr]).T)
	print(df)
	df.columns=['z_inter','min_dist'+title]
	return df
	"""end of comment"""





"""Read RS"""

# main_path=pa.main_path

# plots_folder=pa.plots_folder

# RS_name="RS.BodyCT1_CBCT_1_6_11_16_20_23.dcm"

# RS_path=main_path+RS_name   

# rs = pydicom.read_file(RS_path)


# x_PTV,y_PTV,z_PTV=fc.get_contour_all_slice("PTV_ALL",RS_path)


"""Get the contours from the Parotid and Thyroid"""
# x_p,y_p,z_p=fc.get_contour_all_slice("Parotid",RS_path)
# x_t,y_t,z_t=fc.get_contour_all_slice("Thyroid",RS_path)

# print("min_z_p",min(z_p))

# print("max_z_p",max(z_p))


# print("min_z_t",min(z_t))

# print("max_z_t",max(z_t))

"""Sort contours"""

#x_PTV,y_PTV,z_PTV=sort_contours(x_PTV,y_PTV,z_PTV)#This is not sorted nicely

#body_names=["BODY","Body-1","Body-6","Body-11","Body-16","Body-20","Body-23"]

#body_names=["Body-20"]


# def get_cont_all_fractions(body_names,RS_path,PTV_name,ID):


def get_cont_all_fractions(body_names,RS_path,x_PTV,y_PTV,z_PTV,ID):





	df = pd.DataFrame(columns = ['z_inter','min_dist'])
	n=len(body_names)

	anterior_half_list=[]

	
	

	for i in range(n):



		print("\n\n body name \n",body_names[i])



		x_CT,y_CT,z_CT=fc.get_contour_all_slice(body_names[i],RS_path)

		# x_PTV,y_PTV,z_PTV=fc.get_contour_all_slice(PTV_name,RS_path)


		x_CT,y_CT,z_CT=sort_contours(x_CT,y_CT,z_CT)

		# x_PTV,y_PTV,z_PTV=sort_contours(x_PTV,y_PTV,z_PTV)

		"""Get max contour"""


		x_CT,y_CT,z_CT=max_cnt_CT(x_CT,y_CT,z_CT)


		

		#print(i)

		x0,y0,x1,y1,inter1=all_conts_two_structures(x_CT,y_CT,z_CT,x_PTV,y_PTV,z_PTV)
		#print(inter1)

		# print("53 slice y0",y0[53])

		# print("len ant",len(anterior_half_list))

		# print("anterior half list SLICE 53",anterior_half_list)

		"""The following part is to take the anterior part"""

		# if i==0:

		# 	anterior_half_list=get_ant_half_all_slices(y0)

		# x0,y0=apply_ah_all_slices(x0,y0,anterior_half_list)
		# x1,y1=apply_ah_all_slices(x1,y1,anterior_half_list)

		#print("53 slice",x0[53],y0[53])


		df_CT=min_distance_all_slices(x0,y0,x1,y1,inter1,body_names[i])
		"""I will comment this part"""

		df_CT['z_inter']=df_CT['z_inter'].astype("float")

		#print(df_CT)

		# df_CT['Body']=body_names[i]

		if len(df)==0:
			df=df_CT


			#df.Weight.astype('int64')
		else:



			df=pd.merge(df,df_CT,how='inner', on='z_inter')

		df["ID"]=ID

		print(df)

	return df
	"""End of comment"""


def get_Body_ROI_names(RS_path):

	rs = pydicom.read_file(RS_path)

	body_names=[]
	for key in rs.StructureSetROISequence:
		
		ROI=key.ROIName

		#print(ROI)

		if (ROI[0:4]=="BODY") | (ROI[0:4]=="Body"):

			body_names.append(ROI)

	return list(np.unique(body_names))


# ids_list=pa.ids_list

# patients=list(pd.read_csv(ids_list)["PatientId"])

# pat_id=patients[0]

# path=pa.main_path+str(pat_id)+'/'

# RS_path=path+fc.find_RS_files(path)

# rs = pydicom.read_file(RS_path)

# plots_folder=pa.plots_folder



# body_names=get_Body_ROI_names(RS_path)
# body_names=["BODY"]
# #body_names=["Body-3"]

# #body_names.remove("Body-1")

# print(body_names)

# x_PTV,y_PTV,z_PTV=fc.get_contour_all_slice("PTV_ALL",RS_path)

# x_PTV,y_PTV,z_PTV=sort_contours(x_PTV,y_PTV,z_PTV)

# df=get_cont_all_fractions(body_names,RS_path,x_PTV,y_PTV,z_PTV,str(pat_id))

# df.to_csv("min_dist"+str(pat_id))


