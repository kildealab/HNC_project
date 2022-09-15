from patients_list import *
from helpers import *
from param_xvalues import *

PATH_DEST = '/home/james/PTV-body_distance/07_streamlined-pipeline/Data/line_segments/'
ROWS = ['x_pin', 'y_pin', 'z_pin', 'x_pout', 'y_pout', 'z_pout', 'body_out?']

def get_xmin_points_from_contours(contours_PTV, contours_body, r_frac=RADIUS_FRAC, smooth_iter=0):
    # ================================================================================
    contours_body, contours_PTV = trim_contours_to_match_z(contours_body, contours_PTV)
    cloud_PTV = pv.PolyData(contours_PTV)
    cloud_body = pv.PolyData(contours_body)
    surface_body = get_surface_marching_cubes(contours_body).smooth(n_iter=smooth_iter)
    inner_PTV, outer_PTV = split_cloud_by_surface(cloud_PTV, surface_body)
    xdiff, ydiff = get_bounding_box_dimensions(contours_body)
    if len(outer_PTV.points) > 0 and xdiff > ydiff:
        num_outer_points_ref = len(outer_PTV.points)
        cloud_PTV_trim, h, k, r = trim_posterior_PTV(cloud_PTV, contours_body, 1)
        inner_PTV_new, outer_PTV_new = split_cloud_by_surface(cloud_PTV_trim, surface_body)
        num_outer_points_new = len(outer_PTV_new.points)
        if num_outer_points_new < num_outer_points_ref:
            cloud_PTV_trim, h, k, r = trim_posterior_PTV(cloud_PTV, contours_body, r_frac)
            inner_PTV, outer_PTV = split_cloud_by_surface(cloud_PTV_trim, surface_body)
    # ================================================================================
    
    # get points in the xmin line segment
    if len(outer_PTV.points) > 0:
        idx_PTV, idx_cloud = get_index_farthest_point_from_cloud(outer_PTV, cloud_body)
        pin_xmin = cloud_body.points[idx_cloud]
        pout_xmin = outer_PTV.points[idx_PTV]
        idx_pbody = 0 # determines whether body point is outer or inner point
    else:
        idx_PTV, idx_cloud = get_index_nearest_point_from_cloud(inner_PTV, cloud_body)
        pin_xmin = inner_PTV.points[idx_PTV]
        pout_xmin = cloud_body.points[idx_cloud]
        idx_pbody = 1
        
    return pin_xmin, pout_xmin, idx_pbody


def pipeline_line_segments(param_name='xmin_line-segments'):
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
                    pin_xmin, pout_xmin, idx_pbody = get_xmin_points_from_contours(contours_PTV, contours_body)
                    print('\tProcess time of parameters (' + key_body + '): ' + str(process_time()-t1) + ' s')
                except:
                    print('\tAn error occured while calculating the parameters for patient ' + str_pat_id + ' with body key ' + key_body + '.')
                    print('\tSkipping ' + key_body + '.')
                    continue
                
                # RECORD PARAMETERS
                # keep order same as ROWS!
                params = []
                params.append(pin_xmin[0])
                params.append(pin_xmin[1])
                params.append(pin_xmin[2])
                params.append(pout_xmin[0])
                params.append(pout_xmin[1])
                params.append(pout_xmin[2])
                params.append(idx_pbody)
                # record calculated values under body key
                df[key_body] = params
                # ================================================================================
                
            # write data to csv
            df.to_csv(out_path, index=False)
            print('\t' + param_name + ' printed to csv: ' + out_path)
            print('Elapsed time for patient: ' + str((process_time()-t0)/60) + ' min')
    print('DONE! Elapsed time for pipeline: ' + str((process_time()-t_init)/3600) + ' hours')
    

if __name__ == "__main__":
    pipeline_line_segments()