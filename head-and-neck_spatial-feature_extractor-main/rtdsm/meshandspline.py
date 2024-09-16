import numpy as np
import pyvista as pv
from rtdsm.helpers import find_nearest
# packages required to create structure meshes
from skimage.measure import marching_cubes
from skimage.draw import polygon
from scipy.spatial import distance

def sample_centerline(CenterLine,step,stepType='dist'):
    """
    Given a spline defining a path through an ROI, defines equidistant points
    along the centerline for sampling.
    
    Parameters
    ----------
    CenterLine : pyvista.PolyData
        A densely defined spline of a path through the length of an ROI.
    step : float or int
        Depending on the step type, specifies either the desired distance between
        points or the total number of points to define along the centerline. 
    stepType : str, optional     
        Specifies if sampling should be performed with a set step size or set total
        number of points. Valid entries are 'dist' and 'nsteps'. Defaults to 'dist'.

    Returns
    -------
    sampinds : numpy.ndarray
        Indices of points in the CenterLine.points object at which to sample.  
    """
    if stepType not in ('dist','nsteps'):
        raise Exception('The stepType given is not a valid option.')
    #STEP1: get total length of the centerline spline
    totlen = CenterLine['arc_length'][-1]
        #STEP1b: if using nsteps method, figure out how big the step size should be
    if stepType=='nsteps':
        step = totlen/step
    #STEP2: set the sampling points along the spline and get the sampling points
    sampdists = np.arange(0,totlen,step)
    sampinds = [] 
    for samp in sampdists:
        sampinds.append(find_nearest(CenterLine['arc_length'],samp))
    #STEP3: return the indicies of the sampling points (we give this instead of just the points themseves so can easily use to get points and path lengths at each sampling point)
    return sampinds

def get_spline_tangents(Spline):
    """
    Calculates tangent vectors for each point along a pyvista spline.

    This requires the spline's gradient to have been calculated with the
    compute_derivative() method.

    Parameters
    ----------
    Spline : pyvista.PolyData     
        A pyvista spline objected returned from the pyvista compute_derivative() 
        method.

    Returns
    -------
    tangent :  numpy.ndarray    
        Array of tangent vectors for all points along the spline (shape M x 3).
    """
    #STEP1: calculate the gradient of the spline at all points and normalize to unit vectors
    Spline['gradient'] = Spline['gradient']/(np.linalg.norm(Spline['gradient'],axis=1)[:, np.newaxis] )
    #STEP2: Calculate the vector of a line passing between two points that neighbour each point on the spline
        #obviously skipping the first and last points, which is fine because we fix their sliceing normals as (0,0,1)
        #This will decide the direction of the line in the tangent plane
    StepVects = np.zeros((len(Spline.points),3))
    StepVects[1:-1] = Spline.points[2:] - Spline.points[:-2]
    StepVects[0],StepVects[-1] = [0,0,1],[0,0,1]
    #STEP3: project the vectors onto the tangent plane to get the final spline tangents
    OutOfPlaneComponent = np.sum(StepVects*Spline['gradient'], axis=1)[:, np.newaxis]*Spline['gradient'] #dot product times the unit normal gradient vectors
    tangent = (StepVects - OutOfPlaneComponent)/(np.linalg.norm((StepVects - OutOfPlaneComponent),axis=1)[:, np.newaxis] ) #point to point vector minus the out of plane component, normalized
    tangent[0],tangent[-1] = [0,0,1],[0,0,1] #set the first and last slices to defaults
    return tangent

def get_cubemarch_surface(PointArray,VoxSize):
    """
    Creates a surface mesh of an ROI from a pointcloud dataset using the marching
    cubes algorithm.
                
    Parameters
    ----------
    PointArray : numpy.ndarray     
        An array of point coordinates representing the vertices of the desired mesh.
    VoxSize : numpy.ndarray or list       
        The dimensions (x,y,z) of a voxel in the image the ROI was generated from.
        Required to create the mask used by the matching cube algorithm.

    Returns
    -------
    MmVerts : numpy.ndarray        
        An array of the vertices of the final mesh.
    faces : numpy.ndarray         
        A 2D array specifying the faces of the mesh. Each row lists the indices
        of the vertices that comprise one face.
    pvFaces : numpy.ndarray         
        The same data as faces, just pre-formatted in a pyvista readable array 
        (each row begins with the number of vertices that define the face).

    Notes
    -----
    Meshes created with a marching cubes algorithm will not include holes (like some Delaunay
    meshes do that can cause issues for other functions in this package), but will be 
    somewhat voxelated. For this reason it is recommended to apply smoothing to the ROI mesh
    once formatted as a pyvista.PolyData object. 
    """
    #update 2022-01-19: changed definition of nSlices and Slices to not be sensitive to very small changes in Z -HP
    #STEP1: Get the min X,Y values to set the top corner of the slices
    Xmin,Ymin = np.nanmin(PointArray[:,0]),np.nanmin(PointArray[:,1])
    Zmin = np.nanmin(PointArray[:,2])
    CornerOrg = [Xmin - 2*VoxSize[0], Ymin - 2*VoxSize[1]]
    #STEP2: convert the XY values to index positions in a 3D array 
    PointArray[:,0] = np.round((PointArray[:,0]-CornerOrg[0])/VoxSize[0])
    PointArray[:,1] = np.round((PointArray[:,1]-CornerOrg[1])/VoxSize[1])
    #STEP3: Make a list of the indices where the slice Z value changes by at least half the Z resolution
    deltaslice = PointArray[:,2] - np.roll(PointArray[:,2],1)
    Slices = np.where((np.abs(deltaslice) >= VoxSize[2]/2)|(np.isnan(deltaslice)) ) [0] 
    #STEP4: Determine how many slices are needed to cover the full structure and make an empty grid
    uniqueSlices = np.unique(PointArray[Slices,2][~np.isnan(PointArray[Slices,2])])
    nSlices = len(uniqueSlices) 
    GridMaxInd = np.nanmax(PointArray[:,:2])   #the max of the X and Y index values
    MaskGrid = np.zeros((int(nSlices),int(GridMaxInd +2),int(GridMaxInd+2))) #NOTE: using ZYX here
    #STEP5: Loop through the slices and add the information to the mask
    for i in range(len(Slices)):
        CurrentSlice = PointArray[Slices[i],2]
        if np.isnan(CurrentSlice):
            continue #if a seperator nan entry, skip over
        #get the list of points for that polygon
        end = -1 if i == len(Slices)-1 else Slices[i+1]
        iPoints = PointArray[Slices[i]:end,:2] #just need the X and Y points
        if len(iPoints) == 0:### TODO
            continue
        #Make a polygon mask from the points
        rr, cc = polygon(iPoints[:,1], iPoints[:,0])        #r AKA row is Y, c AKA col is X
        sliceInd = np.where(uniqueSlices == CurrentSlice)[0]
        MaskGrid[sliceInd,rr, cc] = 1 #add it to the 3d mask grid
    #STEP6: using the 3D mask, run the marching cube algorithm
    verts, faces, normals, values = marching_cubes(MaskGrid)
    #STEP7: Convert the vertices back to mm in the patient coordinate system
    MmVerts = np.zeros(np.shape(verts))
    MmVerts[:,2] =  verts[:,0]*VoxSize[2]+ np.nanmin(PointArray[:,2])#convert Z inds to mm positions
    MmVerts[:,1] =  verts[:,1]*VoxSize[1]+ CornerOrg[1] #convert Y inds to mm positions
    MmVerts[:,0] =  verts[:,2]*VoxSize[0]+ CornerOrg[0] #convert X inds to mm positions
    # also, create an array of face information in a pyvista readable format
    pvFaces = np.zeros((len(faces),4),dtype=np.uint32)
    pvFaces[:,0],pvFaces[:,1],pvFaces[:,2],pvFaces[:,3] = 3, faces[:,0], faces[:,1], faces[:,2] 
    return MmVerts, faces, pvFaces


def get_cubemarch_surface_clusters(PointArray,VoxSize, max_away=20.):
    #update 2022-01-19: changed definition of nSlices and Slices to not be sensitive to very small changes in Z -HP
    #STEP1: Get the min X,Y values to set the top corner of the slices
    Xmin,Ymin = np.nanmin(PointArray[:,0]),np.nanmin(PointArray[:,1])
    Zmin = np.nanmin(PointArray[:,2])
    CornerOrg = [Xmin - 2*VoxSize[0], Ymin - 2*VoxSize[1]]
    #STEP2: convert the XY values to index positions in a 3D array 
    PointArray[:,0] = np.round((PointArray[:,0]-CornerOrg[0])/VoxSize[0])
    PointArray[:,1] = np.round((PointArray[:,1]-CornerOrg[1])/VoxSize[1])
    #STEP3: Make a list of the indices where the slice Z value changes by at least half the Z resolution
    deltaslice = PointArray[:,2] - np.roll(PointArray[:,2],1)
    Slices = np.where((np.abs(deltaslice) >= VoxSize[2]/2)|(np.isnan(deltaslice)) ) [0] 
    #STEP4: Determine how many slices are needed to cover the full structure and make an empty grid
    uniqueSlices = np.unique(PointArray[Slices,2][~np.isnan(PointArray[Slices,2])])
    nSlices = len(uniqueSlices) 
    GridMaxInd = np.nanmax(PointArray[:,:2])   #the max of the X and Y index values
    MaskGrid = np.zeros((int(nSlices),int(GridMaxInd +2),int(GridMaxInd+2))) #NOTE: using ZYX here
    #STEP5: Loop through the slices and add the information to the mask
    for i in range(len(Slices)):
        CurrentSlice = PointArray[Slices[i],2]
        sliceInd = np.where(uniqueSlices == CurrentSlice)[0]
        if np.isnan(CurrentSlice):
            continue #if a seperator nan entry, skip over
        #get the list of points for that polygon
        if i == len(Slices)-1:
            iPoints = PointArray[Slices[i]:,:2]
        else:
            iPoints = PointArray[Slices[i]:Slices[i+1],:2] #just need the X and Y points
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
            MaskGrid[sliceInd,rr, cc] = 1 #add it to the 3d mask grid
    #STEP6: using the 3D mask, run the marching cube algorithm
    verts, faces, normals, values = marching_cubes(MaskGrid)
    #STEP7: Convert the vertices back to mm in the patient coordinate system
    MmVerts = np.zeros(np.shape(verts))
    MmVerts[:,2] =  verts[:,0]*VoxSize[2]+ np.nanmin(PointArray[:,2])#convert Z inds to mm positions
    MmVerts[:,1] =  verts[:,1]*VoxSize[1]+ CornerOrg[1] #convert Y inds to mm positions
    MmVerts[:,0] =  verts[:,2]*VoxSize[0]+ CornerOrg[0] #convert X inds to mm positions
    # also, create an array of face information in a pyvista readable format
    pvFaces = np.zeros((len(faces),4),dtype=np.uint32)
    pvFaces[:,0],pvFaces[:,1],pvFaces[:,2],pvFaces[:,3] = 3, faces[:,0], faces[:,1], faces[:,2] 
    return MmVerts, faces, pvFaces