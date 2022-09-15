from helpers import *

from time import process_time
from scipy.stats import sem
from scipy.spatial import KDTree
from skimage.measure import marching_cubes

import gc

PATH_DEST = '/home/james/PTV-body_distance/07_streamlined-pipeline/Data/xvalues/'
ROWS = ['xmin', 'xmax', 'xlen', 'xmed', 'xave', 'xstd', 'xsem']

def get_index_farthest_point_from_cloud(points, cloud):
    tree = KDTree(cloud.points)
    dist_points, idxes_cloud = tree.query(points.points)
    max_d = max(dist_points)
    idx_point = np.where(dist_points == max_d)[0][0]
    idx_cloud = idxes_cloud[idx_point]
    return idx_point, idx_cloud

def get_index_nearest_point_from_cloud(points, cloud):
    tree = KDTree(cloud.points)
    dist_points, idxes_cloud = tree.query(points.points)
    min_d = min(dist_points)
    idx_point = np.where(dist_points == min_d)[0][0]
    idx_cloud = idxes_cloud[idx_point]
    return idx_point, idx_cloud

def get_distances_of_points_from_cloud(points, cloud):
    tree = KDTree(cloud.points)
    dist_points, idxes_cloud = tree.query(points.points)
    return dist_points

def get_distances_from_contours(contours_PTV, contours_body, r_frac=RADIUS_FRAC, smooth_iter=0):
    # ================================================================================
    # trim body and PTV contours to have matching z limits
    contours_body, contours_PTV = trim_contours_to_match_z(contours_body, contours_PTV)
    
    # get pyvista objects of PTV and body
    cloud_PTV = pv.PolyData(contours_PTV)
    cloud_body = pv.PolyData(contours_body)
    
    # get body surface
    surface_body = get_surface_marching_cubes(contours_body).smooth(n_iter=smooth_iter)
    
    # split PTV cloud into inner and outer
    inner_PTV, outer_PTV = split_cloud_by_surface(cloud_PTV, surface_body)
    xdiff, ydiff = get_bounding_box_dimensions(contours_body)
    
    if len(outer_PTV.points) > 0 and xdiff > ydiff:
        # trim PTV cloud that's protruding outside of FOV
        num_outer_points_ref = len(outer_PTV.points)
        cloud_PTV_trim, h, k, r = trim_posterior_PTV(cloud_PTV, contours_body, 1)
        # redefine outer_PTV and inner_PTV
        inner_PTV_new, outer_PTV_new = split_cloud_by_surface(cloud_PTV_trim, surface_body)
        num_outer_points_new = len(outer_PTV_new.points)
        
        # if something was trimmed, trim again using smaller bounding cylinder
        if num_outer_points_new < num_outer_points_ref:
            cloud_PTV_trim, h, k, r = trim_posterior_PTV(cloud_PTV, contours_body, r_frac)
            inner_PTV, outer_PTV = split_cloud_by_surface(cloud_PTV_trim, surface_body)
    # ================================================================================
        
    if len(outer_PTV.points) > 0:
        distances_outer = -1*get_distances_of_points_from_cloud(outer_PTV, cloud_body)
        distances_inner = get_distances_of_points_from_cloud(inner_PTV, cloud_body)
        distances = np.concatenate((distances_outer,distances_inner))
    else:
        distances = get_distances_of_points_from_cloud(inner_PTV, cloud_body)
    return distances

def pipeline_xvalues(param_name='xvalues'):
    existing_patients = [csv_filename.split('_')[-1].split('.')[0] for csv_filename in os.listdir(PATH_DEST)]
    t_init = process_time()
    for str_pat_id in PATIENTS_ALL:
        # check if patient already has csv
        if str_pat_id in existing_patients:
            print('Patient already has csv:' + str_pat_id)
            continue
        else:
            print('Processing patient: ' + str_pat_id)
            t0 = process_time()
            # get path to RS file
            path_RS = get_path_RS(str_pat_id,PATH_SRC)

            # initialize dataframe and define output file name
            df = pd.DataFrame({param_name : ROWS})
            out_file = param_name + '_' + str_pat_id + '.csv'
            out_path = os.path.join(PATH_DEST,out_file)
            
            # get body keys and sort
            keys_body = get_body_keys(path_RS)
            sorted_keys_body = sort_body_keys(keys_body)
            
            # get contours of PTV
            contours_PTV = rtdsm.get_pointcloud('PTV_ALL', path_RS, False)[0]
            for key_body in sorted_keys_body:
                # get contours of body corresponding to body key
                contours_body = rtdsm.get_pointcloud(key_body, path_RS, False)[0]
                
                # skip body key if contour array is empty
                if len(contours_body) == 0:
                    print('\tSkipping ' + key_body + '. Body contours array is empty.')
                    continue
                
                # ================================================================================
                # CALCULATE PARAMETERS
                try:
                    t1 = process_time()
                    distances = get_distances_from_contours(contours_PTV, contours_body)
                    print('\tProcess time of parameters (' + key_body + '): ' + str(process_time()-t1) + ' s')
                except:
                    print('\tAn error occured while calculating the parameters for patient ' + str_pat_id + ' with body key ' + key_body + '.')
                    print('\tSkipping ' + key_body + '.')
                    continue
                    
                # RECORD PARAMETERS
                # keep order same as ROWS!
                params = []
                params.append(np.min(distances))
                params.append(np.max(distances))
                params.append(np.size(distances))
                params.append(np.median(distances))
                params.append(np.mean(distances))
                params.append(np.std(distances))
                params.append(sem(distances))
                # record calculated values under body key
                df[key_body] = params
                # ================================================================================
                
            # write data to csv
            df.to_csv(out_path, index=False)
            print('\t' + param_name + ' printed to csv: ' + out_path)
            print('Elapsed time for patient: ' + str((process_time()-t0)/60) + ' min')
    print('DONE! Elapsed time for pipeline: ' + str((process_time()-t_init)/3600) + ' hours')
    

if __name__ == "__main__":
    pipeline_xvalues()