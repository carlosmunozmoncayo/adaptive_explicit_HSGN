import numpy as np
import os
import sys

try:
    from clawpack.amrclaw import region_tools
except:
    raise Exception("Issue importing region_tools from clawpack.amrclaw.")

def masks_from_clusters(ref_mask,clusters,n_clusters,points):
    """
    This function creates masks from the clusters of points.
    Input:
    ref_mask: array-like
        A 2D array of shape (m,n) containing the reference mask.
    clusters: array-like
        An array of shape (n,) containing the cluster number for each point.
    n_clusters: int
        The number of clusters found.
    points: array-like
        An array of shape (n,2) containing the x and y coordinates of the points.
    Output:
    cluster_masks: list of arrays
        A list of 2D arrays of shape (m,n) containing the masks for each cluster.
    """
    cluster_masks = []
    for i in range(n_clusters):
        cluster_mask = np.zeros_like(ref_mask)
        cluster_mask[points[clusters==i+1][:,0],points[clusters==i+1][:,1]] = True
        #For clusters with a single point, pad the cluster with a 3x3 square around the point
        if cluster_mask.sum() == 1:
            idx = np.argwhere(cluster_mask)
            cluster_mask[idx[0][0]-1:idx[0][0]+2,idx[0][1]-1:idx[0][1]+2] = True
            print(f"Padding cluster {i}")
        cluster_masks.append(cluster_mask)
    return cluster_masks

def ruled_rects_from_masks(clustered_masks, topo, ixy="y", method=0,
                           list_ixy=None, list_methods=None):
    """
    This function creates ruled rectangles from the clustered masks.
    Input:
    clustered_masks: list of arrays
        A list of 2D arrays of shape (m,n) containing the masks for each cluster.
    topo: Topography object
        A Topography object from clawpack.amrclaw.
    ixy: string (default "y")
        The direction of the ruled rectangle.
        "x" for horizontal ruled rectangles.
        "y" for vertical ruled rectangles.
    method: int (default 0)
    list_ixy: list of strings (default None)
        A list of strings with the direction of the ruled rectangle for each cluster.
    list_methods: list of ints (default None)
        A list of ints with the method for the ruled rectangle for each cluster.
    Output:
    list_ruled_rects: list of RuledRectangle objects
    """
    if list_ixy is None:
        list_ixy = [ixy]*len(clustered_masks)
    if list_methods is None:
        list_methods = [method]*len(clustered_masks)
    list_ruled_rects = []
    for i in range(len(clustered_masks)):
        #Get direction of rr
        cluster_mask = clustered_masks[i]
        #Padding for single point clusters
        if cluster_mask.sum() == 1:
            idx = np.argwhere(cluster_mask)
            cluster_mask[idx[0][0]-1:idx[0][0]+2,idx[0][1]-1:idx[0][1]+2] = True
            print(f"Padding cluster {i}")
        ixy = list_ixy[i]
        method = list_methods[i]
        #Create the ruled rectangle
        rrect = region_tools.ruledrectangle_covering_selected_points(topo.X, topo.Y,
                                                                    cluster_mask, ixy=ixy, method=method,
                                                                    verbose=False)
        list_ruled_rects.append(rrect)
    return list_ruled_rects

def improve_convexity(rr,max_iter=100):
    iter = 0
    improved = True
    size_rr = len(rr.s)
    while iter<max_iter and improved:
        iter+=1
        #Replace nonconvex local extrema by weighted averages of neighbors
        left_upper = rr.upper[:-2]; right_upper = rr.upper[2:]
        left_lower = rr.lower[:-2]; right_lower = rr.lower[2:]
        left_s = rr.s[:-2]; right_s = rr.s[2:]
        local_minima = np.logical_and(rr.upper[1:-1]<=left_upper,rr.upper[1:-1]<=right_upper)
        local_maxima = np.logical_and(rr.lower[1:-1]>=left_lower,rr.lower[1:-1]>=right_lower)
        linear_interp_upper = (right_s-rr.s[1:-1])/(right_s-left_s)*left_upper+(rr.s[1:-1]-left_s)/(right_s-left_s)*right_upper
        linear_interp_lower = (right_s-rr.s[1:-1])/(right_s-left_s)*left_lower+(rr.s[1:-1]-left_s)/(right_s-left_s)*right_lower
        rr.upper[1:-1] = np.where(local_minima, linear_interp_upper, rr.upper[1:-1])
        rr.lower[1:-1] = np.where(local_maxima, linear_interp_lower, rr.lower[1:-1])
        improved = local_minima.any() or local_maxima.any()


if __name__ == "__main__":
    print("Nothing to be done from more_region_tools.py")