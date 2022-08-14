"""
Splicing Forgeries -- 
This module implements a splicing

-------------------------------------
Author: Phillipe Cardenuto
email: phillipe.cardenuto@gmail.com
--------------------------------------
Version:                           1.1
Version Date:               Aug , 2020

This is file is part of Forgery Library
"""


import cv2
from ..utils import _get_start_x_y_from_slice, _find_background
from ..utils import _select_random_objs 
from .copy_move import _border_forgery_post_processing, _copy_move
from scipy import ndimage as ndi
import numpy as np
from PIL import Image
from collections import defaultdict



def splicing_forgery(donor, donor_map, host, host_map, n_objects=1):
    """
    Apply a splicing forgery from donor to host using n_objecrs form donor.
    The method selects an object from the donor using the donor_map 
    and paste it on a background region of host. To find a suitable background region
    the method uses the host_map.
    If it is not possible to select a background region that any object from the donor fits
    ,then then the returning maps will containg only zeros.
        
    Warning
    -------
    Make sure to pass a uint8 image
    
    Parameters
    ----------
    donor : image opencv obj
        Source of the object that will be placed on im_dest
    donor_map: A binary image in opencv format
        This binary map indicates where the objects are in the im_src
        255 -> object
        0   -> background
    host: image opencv obj
        Image that will be tampered
    host_map: A binary image in opencv format
        This binary map indicates where the objects are in the im_src
        255 -> object
        0   -> background
    n_objects: int (optional) = 1
        Nuber of objects envolved on the splicing
        
    Return
    ------
    f_img: opencv Image
        Splicing result image
    donor_splicing_map: opencv Image
        A binary map that indicate the region that the splicing occurs on donor
    host_splicing_map: opencv Image
        A binary map that indicate the region that the splicing occurs on host
    """
    # Select n_objects from the donor_map 
    selected_donor_objs_mask, _ = _select_random_objs(donor_map.copy(),
                                                 n_objects,
                                                 remove_obj_border=True)
    
    # Dilate the selected objects to cover all its borders
    selected_donor_objs_mask[selected_donor_objs_mask > 0 ] = 255
    kernel = np.ones((7,7),np.uint8)
    selected_donor_objs_mask = cv2.dilate(selected_donor_objs_mask,\
                                          kernel,iterations = 1)
    # Label each selected donor object
    selected_donor_objs_mask, _ = ndi.label(selected_donor_objs_mask)
    objs_bbox =  ndi.find_objects(selected_donor_objs_mask)
    
    # Get a copy of the host image in which the forgery will
    # be applied
    f_img = host.copy()
    # Initialize groundtruth map and gt object index
    gt_index = 1 
    donor_gt_map = np.zeros_like(donor_map).astype(np.uint8)
    host_gt_map = np.zeros_like(host_map).astype(np.uint8)
    
    # Loop into each bbox from the donor
    # trying to find a region on host to fit it.
    for obj_id, obj_loc in enumerate(objs_bbox,start=1):
        
        # Bounding Box height and width
        bb_h, bb_w = selected_donor_objs_mask[obj_loc].shape[0:2]
        
        # Find a region that the object could fit
        # This region can not overlap any already selected region
        # or any object
        forbidden_area = cv2.bitwise_or(host_map,host_gt_map)
        forbidden_area[forbidden_area>0] = 255
        background_region = _find_background(forbidden_area,(bb_h,bb_w))
    
        # If there is none background region, the selected object
        # doesn't fit in any background. It is not possible to paste
        # this object into host. 
        if background_region is None:
            continue
    
        # Get the region of the candidate object
        omask = np.isin(selected_donor_objs_mask, [obj_id])
        omask[omask>0] = 255
        omask = omask.astype(np.uint8)
        
        # Apply the splicing (copy_move) forgery
        f_img, donor_gt_aux, host_gt_aux = _copy_move(donor,omask,f_img,background_region)
        
        # Update host and donor groundtruth
        host_gt_aux[host_gt_aux>0] = gt_index 
        host_gt_map = cv2.bitwise_or(host_gt_aux,host_gt_map)
        donor_gt_aux[donor_gt_aux>0] = gt_index 
        donor_gt_map= cv2.bitwise_or(donor_gt_aux,donor_gt_map)
        
        # update gt_index
        gt_index +=1
    
    # Apply post-processing
    f_img = _border_forgery_post_processing(f_img, host_gt_map)
    
    return f_img, donor_gt_map, host_gt_map
