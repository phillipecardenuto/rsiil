"""
Retouching Forgeries -- Cleaning
This module implements cleaning with two different approaches
- Erasing the object using a selected background
- Erasing the object using a blending background
- Erasing the object using a inpainting
    Exemplar-Based Inpainting Criminisi et al (2003)
    

-------------------------------------
Author: Phillipe Cardenuto
email: phillipe.cardenuto@ic.unicamp.br
--------------------------------------
Version:                           1.1
Version Date:               Sep , 2020

This is file is part of Forgery Library
"""

import cv2
from ..utils import _select_random_objs, _find_background
from ..utils import _get_start_x_y_from_slice
from ..duplication.copy_move import _copy_move
from ..retouching import retouching_blurring
from .inpaint import Inpainter
from scipy import ndimage as ndi
import numpy as np
import random
from PIL import Image


def _bbox_enlarge(image_shape, obj_loc, enlarge_size = 100):
    """
    Enlarge a bounding box wihtout falling outside an image
    """
    height, width = image_shape[:2]
    crop_width = min(enlarge_size+ obj_loc[1].stop - obj_loc[1].start, width)
    crop_height = min(enlarge_size+ obj_loc[0].stop - obj_loc[0].start, height)

    x0 = max(0, obj_loc[1].start - enlarge_size)
    y0 = max(0, obj_loc[0].start - enlarge_size)

    x1 = x0 + crop_width
    y1 = y0 + crop_height

    return (slice(y0,y1,None), slice(x0,x1,None))


def cleaning_with_inpainting(img, objs_map, n_objects=1, enlarge_size=50 , message=False):
    """
    Remove n_object from 'img', using the inpainting method 
    Object Removal by Exemplar-Based Inpainting Criminisi et al. (2003)
    For a more accurate result, the algorithm selects a rectangle region around
    the mask that cover all the mask.

    This method return the forgery image  along with a groundtruth that 
    indicates which object was removed, and which region was choose to hide 
    the object.
    
    Warning
    -------
    This method might take some time per object to perform inpaiting 

    Parameters
    ----------
    img: Image object on opencv format
    
    objs_map:Image object on opencv format 
        Binary map, indicating the region of all objects on 'img'
    
    n_objects: int
        Maximum number of objects that will be removed.
        The objects are randomly selected.
        This method does not guarantee that n_objects are removed,
        because there may not have enough background area to do this.
        
    enlarge_size: int (optional)
        Enlarge the affected impaiting bounding box area by the input
        If applied to western blots, you might not want to enlarge the area of the impaiting
        by a large number (perhaps 10 is enought).
        
    message: bool (optional)
    Print status of the inpainting process
    
    Return
    ------
    f_img: opencv Image object
        Image with the cleanned applied
    final_gt: opencv Image object
       Image mask with the location of the objects cleanned 
       as its ID
    """


    # Make copy of input image
    f_img = img.copy()
    if len(f_img.shape) >2:
        if f_img.shape[-1] > 3:
            f_img = f_img[:,:,:3]
            #raise Exception("Image must have at most 3 bands.")

    # Select inpainting object candidate
    # Make a working copy of the object
    objs = objs_map.copy()
    # Dilates the working objects to join objects that are too close
    kernel = np.ones((11,11),np.uint8)
    objs = cv2.dilate(objs,kernel,iterations = 2)
    objs[objs>0] = 1

    # Locate joinned objects
    objs, _ = ndi.label(objs)
    obj_location =  ndi.find_objects(objs)

    # Search for candidate objects, i.e objects could be isolated on a rectangle
    # of size (bbx.h +50, bbx.w +50)
    candidate_location = []
    for bbx in obj_location:
        bbx = _bbox_enlarge(objs.shape, bbx, enlarge_size=enlarge_size)
        if len(np.unique(objs[bbx])) == 2:
            candidate_location.append(bbx)

    # If none candidates were found
    if candidate_location == []:
        return img, np.zeros_like(img)

    # Initialize final ground-truth 
    final_gt = np.zeros_like(objs_map)
    selected_bbxs = random.sample(candidate_location, k=min(len(candidate_location), n_objects))

    # Perform cleanning for all selected objects
    for obj_id, obj_loc in enumerate(selected_bbxs,start=1):
        
        #Intialize removal_map
        removal_map = np.zeros_like(objs_map)
        removal_map[obj_loc] = 255
        crop_slice = _bbox_enlarge(img.shape, obj_loc, enlarge_size=enlarge_size)
        
        # Perform inpainting in removal area
        if message:
            print("Processing object -- %d"%obj_id)
        inpainting_res = Inpainter(f_img[crop_slice],
                                   removal_map[crop_slice],
                                   message,
                                   patch_size=25).inpaint()

        # Update image with inpaiting result
        if len(f_img.shape)  == 2:
            inpainting_res = inpainting_res[:,:,0]
        f_img[crop_slice] = inpainting_res

        # Update final ground-truth map
        final_gt[obj_loc] = obj_id
    
    return f_img, final_gt 



def cleaning_using_bf(img, objs_map, n_objects=1, blur_cleaning=True):
    """
    Cleaning using brute force
    Remove n_object from 'img', pasting on top of each a contiguos random region from the background. After this
    a Gaussian filter is applied to obfuscated the noise of the copied region (just like a human would do).
    During the multiple removals, the same background region
    is not used twice.
    
    This method return the forgery image  along with a groundtruth that indicates which object was removed, and
    which region was choose to hide the object.
    
    Warning
    -------
    Make sure to pass a uint8 image
    
    Parameters
    ----------
    img: Image object on opencv format
    
    objs_map:Image object on opencv format 
        Binary map, indicating the region of all objects on 'img'
    
    n_objects: int
        Maximum number of objects that will be removed.
        The objects are randomly selected.
        This method does not guarantee that n_objects are removed,
        because there may not have enough background area to do this.
        
    blur_cleaning: bool
        Apply blurring in the top of cleaned regions
    
    Return
    ------
    f_img: opencv Image object
        Image with the cleanned applied
    bg_gt : opencv Image object
        Mask of the selected background region to hide the obj which the id is ID
    obj_gt: opencv Image object
        Maks of each object cleanned from the image,
        All objects ID are equal to the ID of the region that cover its
    """
    # Select n_objects from the objs_map  
    selected_objs_mask, _ = _select_random_objs(objs_map.copy(),
                                                 n_objects,
                                                 remove_obj_border=False)
    
    # Dilate the selected objects to cover all its borders
    selected_objs_mask[selected_objs_mask > 0 ] = 255
    kernel = np.ones((11,11),np.uint8)
    selected_objs_mask = cv2.dilate(selected_objs_mask,kernel,iterations = 1)
    # Label each object
    selected_objs_mask, _ = ndi.label(selected_objs_mask)
    objs_bbox =  ndi.find_objects(selected_objs_mask)
    
    # Get Copies from the original inputs
    aux_img = img.copy()
    # Initialize groundtruth map and gt object index
    gt_index = 1 
    bkg_gt_map = np.zeros_like(objs_map).astype(np.uint8)
    obj_gt_map = np.zeros_like(objs_map).astype(np.uint8)
    
    # Loop into each bbox trying to find a region to fit
    # on the respective object
    for obj_id, obj_loc in enumerate(objs_bbox,start=1):
        
        # Bounding Box height and width
        bb_h, bb_w = selected_objs_mask[obj_loc].shape[0:2]
        
        # Find a region that the object could fit
        # This region can not overlap any already selected region
        # or any object
        forbidden_area = cv2.bitwise_or(objs_map,bkg_gt_map)
        background_region = _find_background(forbidden_area,(bb_h,bb_w))

        # If there is none background region, the selected object
        # doesn't fit in any background. It is not possible to clean
        # this object with a contiguos background region.
        if background_region is None:
            continue
    
        # Get the region of the candidate object
        omask = np.isin(selected_objs_mask, [obj_id])
        omask[omask>0] = 255
        omask = omask.astype(np.uint8)
        
        
        # Apply the forgery cleaning, copy the background into the selected object
        aux_img, bkg_gt_aux, obj_gt_aux = _copy_move(img,background_region, aux_img, omask)
        
        if blur_cleaning:
            # Apply blurring on the cleanned area, to smooth the region
            aux_img, _= retouching_blurring(aux_img, obj_gt_aux.copy(),rand=False,gaussian=True)
        
        # Update obj_gt_map, bkg_gt_map, and obj_map
        bkg_gt_aux[bkg_gt_aux>0] = gt_index 
        bkg_gt_map= cv2.bitwise_or(bkg_gt_aux,bkg_gt_map)
        obj_gt_aux[obj_gt_aux>0] = gt_index 
        obj_gt_map= cv2.bitwise_or(obj_gt_aux,obj_gt_map)
        objs_map = cv2.bitwise_or(objs_map,obj_gt_map)
        
        # update gt_index
        gt_index +=1
    
    f_img = aux_img
    return f_img, bkg_gt_map, obj_gt_map 

