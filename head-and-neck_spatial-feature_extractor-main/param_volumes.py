from patients_list import *
from helpers import *

from time import process_time
from skimage.draw import polygon
from scipy.spatial import distance

PATH_DEST = '/home/james/PTV-body_distance/07_streamlined-pipeline/Data/volumes/'
ROWS = ['vol_body', 'vol_PTV_inner', 'vol_PTV_outer']
BIG_AWAY = 400
SMALL_AWAY = 30

# snippet from rtdsm's get_cubemarch_surface
def get_mask_grid(contours, max_away=400.):
    img_res = [IMG_RES[0], IMG_RES[1], get_contour_z_spacing(contours)]
    #STEP1: Get the min X,Y values to set the top corner of the slices
    Xmin,Ymin = np.nanmin(contours[:,0]),np.nanmin(contours[:,1])
    Zmin = np.nanmin(contours[:,2])
    CornerOrg = [Xmin - 2*img_res[0], Ymin - 2*img_res[1]]

    #STEP2: convert the XY values to index positions in a 3D array 
    contours[:,0] = np.round((contours[:,0]-CornerOrg[0])/img_res[0])
    contours[:,1] = np.round((contours[:,1]-CornerOrg[1])/img_res[1])

    #STEP3: Determine how many slices are needed to cover the full structure and make an empty grid
    uniqueSlices = np.unique(contours[:,2][~np.isnan(contours[:,2])])
    nSlices = len(uniqueSlices)
    GridMaxInd = np.nanmax(contours[:,:2])   #the max of the X and Y index values
    MaskGrid = np.zeros((int(nSlices),int(GridMaxInd +2),int(GridMaxInd+2))) #NOTE: using ZYX here
    
    #STEP4: Make a list of the indices where the slice number changes
    deltaslice = contours[:,2] - np.roll(contours[:,2],1)
    Slices = np.where(deltaslice != 0)[0] #indexes where a new polygon begins
    for i in range(len(Slices)):
        CurrentSlice = contours[Slices[i],2]
        sliceInd = np.where(uniqueSlices == CurrentSlice)[0]
        if np.isnan(CurrentSlice):
            continue
        #get the list of points for that polygon
        if i == len(Slices)-1:
            iPoints = contours[Slices[i]:,:2]
        else:
            iPoints = contours[Slices[i]:Slices[i+1],:2] #just need the X and Y points
        # split iPoints into cluster
        idx_breaks = [0]
        if len(iPoints) > 1:
            for i in range(len(iPoints)-1):
                if distance.euclidean(iPoints[i], iPoints[i+1]) < max_away:
                    continue
                else:
                    idx_breaks.append(i+1)
        for i in range(len(idx_breaks)): #Make a polygon mask from the points
            start = idx_breaks[i]
            if i == len(idx_breaks)-1:
                rr, cc = polygon(iPoints[start:,1], iPoints[start:,0])
            else:
                end = idx_breaks[i+1]
                rr, cc = polygon(iPoints[start:end,1], iPoints[start:end,0]) #r AKA row is Y, c AKA col is X
            MaskGrid[sliceInd,rr, cc] = 1
    return MaskGrid

def get_volume(contours, img_res, max_away=400.):
    volume_voxel = img_res[0] * img_res[1] * img_res[2]
    grid_mask = get_mask_grid(contours.copy(), max_away)
    volume = volume_voxel * np.sum(grid_mask)
    return volume

def get_volumes_from_contours(contours_PTV, contours_body, r_frac=RADIUS_FRAC, smooth_iter=0):
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
            cloud_PTV = cloud_PTV_trim
        
    img_res_body = [IMG_RES[0], IMG_RES[1], get_contour_z_spacing(contours_body)]
    img_res_PTV = [IMG_RES[0], IMG_RES[1], get_contour_z_spacing(contours_PTV)]
    
    vol_body = get_volume(contours_body, img_res_body, BIG_AWAY)
    vol_PTV = get_volume(cloud_PTV.points, img_res_PTV, BIG_AWAY)
    if len(outer_PTV.points) > 0:
        vol_PTV_outer = get_volume(outer_PTV.points, img_res_PTV, SMALL_AWAY)
        vol_PTV_inner = vol_PTV - vol_PTV_outer
    else:
        vol_PTV_outer = 0
        vol_PTV_inner = vol_PTV
    return vol_body, vol_PTV_inner, vol_PTV_outer
    
    
def pipeline_volumes(param_name='volumes'):
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
                # try:
                t1 = process_time()
                vol_body, vol_PTV_inner, vol_PTV_outer = get_volumes_from_contours(contours_PTV, contours_body)
                print('\tProcess time of parameters (' + key_body + '): ' + str(process_time()-t1) + ' s')
                # except:
                #     print('\tAn error occured while calculating the parameters for patient ' + str_pat_id + ' with body key ' + key_body + '.')
                #     print('\tSkipping ' + key_body + '.')
                #     continue
                    
                # RECORD PARAMETERS
                # keep order same as ROWS!
                params = []
                params.append(vol_body)
                params.append(vol_PTV_inner)
                params.append(vol_PTV_outer)
                # record calculated values under body key
                df[key_body] = params
                # ================================================================================
                
            # write data to csv
            df.to_csv(out_path, index=False)
            print('\t' + param_name + ' printed to csv: ' + out_path)
            print('Elapsed time for patient: ' + str((process_time()-t0)/60) + ' min')
    print('DONE! Elapsed time for pipeline: ' + str((process_time()-t_init)/3600) + ' hours')
    

if __name__ == "__main__":
    pipeline_volumes()