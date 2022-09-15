import numpy as np
import os
from pyvista import Cylinder

import pandas as pd
import pyvista as pv
import pydicom, rtdsm

from patients_list import *

def get_path_RS(pat_id, path_src):
    path_patient = os.path.join(path_src, pat_id)
    file_RS = [x for x in os.listdir(path_patient) if 'RS' in x][0]
    return os.path.join(path_patient, file_RS)

def get_ROI_keys(RS_file_path):
    RS_file = pydicom.read_file(RS_file_path)
    contour_keys = RS_file.StructureSetROISequence
    return [str(x.ROIName) for x in contour_keys]

def get_body_keys(RS_file_path):
    ROI_keys = get_ROI_keys(RS_file_path)
    return [x for x in ROI_keys if 'body' in x.lower()]

def get_PTV_keys(RS_file_path):
    ROI_keys = get_ROI_keys(RS_file_path)
    return [x for x in ROI_keys if 'ptv' in x.lower()]

def sort_body_keys(keys_body):
    new_keys_body = []
    nums = []
    for key in set(keys_body):
        str_frac_num = key.split('-')[-1]
        if not str_frac_num.lower() == 'body':
            nums.append(int(str_frac_num))
        else:
            new_keys_body.append(key)
    nums = sorted(nums)
    for num in nums:
        for key in keys_body:
            if str(num) == key.split('-')[-1]:
                new_keys_body.append(key)        
    return new_keys_body

def get_contour_z_spacing(contours):
    z_vals = np.array(list(set(contours[:,2])))
    z_vals = z_vals[~(np.isnan(z_vals))]
    sorted_z = np.array(sorted(z_vals))
    diff_arr = sorted_z[:-1] - sorted_z[1:]
    return abs(np.mean(diff_arr))

def trim_contours_to_match_z(contours_1, contours_2): # 1: body, 2: PTV
    spacing_z = get_contour_z_spacing(contours_1)
    max_z = max(contours_1[:,2]) - 3*spacing_z
    min_z = min(contours_1[:,2]) + 3*spacing_z
    contours_1 = np.array([x for x in contours_1 if x[2] < max_z and x[2] > min_z])
    contours_2 = np.array([x for x in contours_2 if x[2] < max_z and x[2] > min_z])
    return contours_1, contours_2

def get_surface_marching_cubes(contours):
    img_res = [IMG_RES[0], IMG_RES[1], get_contour_z_spacing(contours)]
    verts, faces, pvfaces = rtdsm.get_cubemarch_surface(contours.copy(), img_res)
    mesh = pv.PolyData(verts, faces=pvfaces)
    return mesh.extract_surface()

def get_surface_marching_cubes_clusters(contours):
    img_res = [IMG_RES[0], IMG_RES[1], get_contour_z_spacing(contours)]
    verts, faces, pvfaces = rtdsm.get_cubemarch_surface_clusters(contours.copy(), img_res)
    mesh = pv.PolyData(verts, faces=pvfaces)
    return mesh.extract_surface()

def split_cloud_by_surface(cloud, surface):
    cloud.compute_implicit_distance(surface, inplace=True)
    inner_cloud = cloud.threshold(0.0, scalars="implicit_distance", invert=True)
    outer_cloud = cloud.threshold(0.0, scalars="implicit_distance", invert=False)
    return inner_cloud, outer_cloud

# def get_proper_inner_and_outer_PTV(contours_PTV, contours_body, r_frac=RADIUS_FRAC, smooth_iter=0):
#     # trim body and PTV contours to have matching z limits
#     contours_body, contours_PTV = trim_contours_to_match_z(contours_body, contours_PTV)
    
#     # get pyvista objects of PTV and body
#     cloud_PTV = pv.PolyData(contours_PTV)
#     cloud_body = pv.PolyData(contours_body)
    
#     # get body surface
#     surface_body = get_surface_marching_cubes(contours_body).smooth(n_iter=smooth_iter)
    
#     # split PTV cloud into inner and outer
#     inner_PTV, outer_PTV = split_cloud_by_surface(cloud_PTV, surface_body)
#     xdiff, ydiff = get_bounding_box_dimensions(contours_body)
    
#     if len(outer_PTV.points) > 0 and xdiff > ydiff:
#         # trim PTV cloud that's protruding outside of FOV
#         num_outer_points_ref = len(outer_PTV.points)
#         cloud_PTV_trim, h, k, r = trim_posterior_PTV(cloud_PTV, contours_body, 1)
#         # redefine outer_PTV and inner_PTV
#         inner_PTV_new, outer_PTV_new = split_cloud_by_surface(cloud_PTV_trim, surface_body)
#         num_outer_points_new = len(outer_PTV_new.points)
        
#         # if something was trimmed, trim again using smaller bounding cylinder
#         if num_outer_points_new < num_outer_points_ref:
#             cloud_PTV_trim, h, k, r = trim_posterior_PTV(cloud_PTV, contours_body, r_frac)
#             inner_PTV, outer_PTV = split_cloud_by_surface(cloud_PTV_trim, surface_body)
#             cloud_PTV = cloud_PTV_trim
        
#     return inner_PTV, outer_PTV, cloud_PTV

# For a given set of contours and x value, return closest x values that is in one of the contour points
def get_closest_x(x, contours):
    closest_x = contours[0][0]
    min_abs_diff = 10000
    for point in contours:
        current_x = point[0]
        abs_diff = abs(current_x - x)
        if abs_diff < min_abs_diff:
            min_abs_diff = abs_diff
            closest_x = current_x
    return closest_x

# For a given set of contours and x value, return the maximum y value
def get_max_y(x, contours):
    target_x = get_closest_x(x, contours)
    max_y = -1
    for point in contours:
        current_x, current_y = point[0:2]
        if round(current_x,1) == round(target_x,1):
            if current_y > max_y:
                max_y = current_y
    return max_y

# For a given set of contours and x value, return point with maximum y value
def get_point_with_max_y_around_given_x(x, contours):
    target_x = x
    max_y = -1
    for point in contours:
        current_x, current_y = point[0:2]
        if abs(current_x - x) < 2:
            if current_y > max_y:
                max_y = current_y
                target_x = current_x
    return (target_x, max_y)

# For 3 given points in a circle, return the center (h,k) and radius r
def get_h_k_r(point1, point2, point3):
    x1, y1 = point1
    x2, y2 = point2
    x3, y3 = point3
    
    x12 = x1 - x2;
    x13 = x1 - x3;
    y12 = y1 - y2;
    y13 = y1 - y3;
    y31 = y3 - y1;
    y21 = y2 - y1;
    x31 = x3 - x1;
    x21 = x2 - x1;
 
    # x1^2 - x3^2
    sx13 = pow(x1, 2) - pow(x3, 2);
 
    # y1^2 - y3^2
    sy13 = pow(y1, 2) - pow(y3, 2);
 
    sx21 = pow(x2, 2) - pow(x1, 2);
    sy21 = pow(y2, 2) - pow(y1, 2);
 
    f = (((sx13) * (x12) + (sy13) *
          (x12) + (sx21) * (x13) +
          (sy21) * (x13)) // (2 *
          ((y31) * (x12) - (y21) * (x13))));
             
    g = (((sx13) * (y12) + (sy13) * (y12) +
          (sx21) * (y13) + (sy21) * (y13)) //
          (2 * ((x31) * (y12) - (x21) * (y13))));
 
    c = (-pow(x1, 2) - pow(y1, 2) -
         2 * g * x1 - 2 * f * y1);
 
    # eqn of circle be x^2 + y^2 + 2*g*x + 2*f*y + c = 0
    # where centre is (h = -g, k = -f) and
    # radius r as r^2 = h^2 + k^2 - c
    h = -g;
    k = -f;
    sqr_of_r = h * h + k * k - c;
 
    # r is the radius
    r = round(np.sqrt(sqr_of_r), 5);
    return [h, k, r]

def get_bounding_box_dimensions(contours):
    max_x = max(contours[:,0])
    min_x = min(contours[:,0])
    diff_x = max_x - min_x

    max_y = max(contours[:,1])
    min_y = min(contours[:,1])
    diff_y = max_y - min_y
    
    return [diff_x, diff_y]

def trim_posterior_PTV(cloud_PTV, contours_body, r_frac=1):    
    max_x = max(contours_body[:,0])
    min_x = min(contours_body[:,0])

    point0 = (min_x, get_max_y(min_x, contours_body))
    point1 = get_point_with_max_y_around_given_x(min_x/2, contours_body)
    point2 = get_point_with_max_y_around_given_x(0, contours_body)
    point3 = get_point_with_max_y_around_given_x(max_x/2, contours_body)
    point4 = (max_x, get_max_y(max_x, contours_body))

    h1,k1,r1 = get_h_k_r(point0, point1, point4)
    h2,k2,r2 = get_h_k_r(point0, point3, point4)
    # h3,k3,r3 = get_h_k_r(point0, point2, point4)
    h = np.mean([h1,h2])
    k = np.mean([k1,k2])
    r = np.mean([r1,r2])
    
    max_z = max(contours_body[:,2])
    min_z = min(contours_body[:,2])
    z = np.mean([min_z,max_z])
    spacing_z = get_contour_z_spacing(contours_body)
    height = (max_z - min_z) + 2*spacing_z

    bounding_cylinder = Cylinder(center=[h,k,z], direction=[0,0,1], radius=r*r_frac, height=height)
    cloud_PTV.compute_implicit_distance(bounding_cylinder, inplace=True)
    cloud_PTV_trim = cloud_PTV.threshold(0.0, scalars="implicit_distance", invert=True)
    
    return cloud_PTV_trim, h, k, r*r_frac