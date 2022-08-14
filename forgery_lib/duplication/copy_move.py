"""
Copy-Move Forgeries -- 
This module implements a function of copy-move

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
from scipy import ndimage as ndi
import numpy as np
from PIL import Image
from collections import defaultdict


# Auxiliar function
def _copy_move(donor, donor_map, host, host_map):
    """
    Perform a copy move forgery copying the donor region (highlighted on donor_map)
    to the host region (highlighted on host_map). 
    
    The copy move will be applied on just one contiguos region
   
    """
    
    assert donor.shape[:2] == donor_map.shape[:2] , "Donor image must be the same size of Donor Map"
    assert host.shape[:2] == host_map.shape[:2] , "Host image must be the same size of Host Map"
    
    # Resize the Donor image to have the same size of the Host
    orig_donor_shape = donor.shape[:2]
    rows, cols = host.shape[:2]
    donor = cv2.resize(donor,(cols,rows))
    donor_map = cv2.resize(donor_map,(cols,rows))
    
    # Get donor region location
    donor_map[donor_map>0] = 1
    donor_loc = ndi.find_objects(donor_map)[0]
    donor_bb_y, donor_bb_x = _get_start_x_y_from_slice(donor_loc)
    
    # Get Host region location
    host_map[host_map>0] = 1
    host_loc = ndi.find_objects(host_map)[0]
    host_bb_y, host_bb_x = _get_start_x_y_from_slice(host_loc)
    
    # Translate the donor to the position of the host 
    x_shift = host_bb_x - donor_bb_x 
    y_shift = host_bb_y - donor_bb_y 
    
    # Translation Matrix
    M = np.float32([[1, 0, x_shift], [0, 1, y_shift]]) 
    # Get shape of host_map
    rows, cols = host_map.shape[0:2]
    # Donor Translation 
    donor_map_translated = cv2.warpAffine(donor_map.copy(), \
                                  M, dsize=(cols, rows)) 
    
    # Host could be on the edge of the image, so the donor region
    # could appear cropped
    # Assert the shape of the image on the host location
    copy_move_map = cv2.bitwise_and( donor_map_translated, host_map)
    
    # Make the inverse operation to get the real map of  the
    # donor groundtruth region
    M_1 = np.float32([[1, 0, -x_shift], [0, 1, -y_shift]]) 
    donor_gt = cv2.warpAffine(copy_move_map.copy(), M_1, dsize=(cols, rows))
    
    #--------------------------
    # Apply copy_move operation
    #--------------------------
    # Get the Translated donor image
    donor_tranlsated = cv2.warpAffine(donor.copy(), M, dsize=(cols,rows)) 
    # Get a PIL version of the translated donor 
    donor_tranlsated = Image.fromarray(donor_tranlsated)
    
    # Get the copy_move_map on Pil format 
    copy_move_map[copy_move_map>0] = 255
    copy_move_map = Image.fromarray(copy_move_map)
    
    # Get the host image in Pil format
    host_pil = Image.fromarray(host)
        
    # Make Forgery
    host_pil.paste(donor_tranlsated, (0, 0),copy_move_map)
    
    # Transforming host_pil after forgery to numpy object
    f_img = np.array(host_pil)
    
    # Resize Donor Gt to original size
    rows, cols = orig_donor_shape
    donor_gt = cv2.resize(donor_gt,(cols,rows))
    
    # Forgery image, Donor's Groundtruth, Host's Groundtruth
    return f_img, donor_gt, np.array(copy_move_map)


def random_copy_move(img, objs_map, n_objects=1, poisson=False):
    """
    Performs a copy-move selecting at random n_objects from the objs_map
    and pasting it on background regions
    
    Warning
    -------
    Make sure to pass a uint8 image
    
    Parameters
    ----------
    img: Image object on opencv format
    
    objs_map:Image object on opencv format 
        Binary map, indicating the region of all objects on 'img'
    
    n_objects: int
        Maximum number of objects that will be copy-moved.
        The objects are randomly selected.
        This method does not guarantee that n_objects are copy-moved,
        because there may not have enough background area to do this.
    poisson: <bool> optional = False
        Apply poisson blend during the copy-move
    
    Return
    ------
    f_img: opencv Image object
        Image with the cleanned applied
    bg_gt : opencv Image object
        Mask of the selected background region that a objects were pasted.
    obj_gt: opencv Image object
        Maks of each object copy-moved from the image,
        All objects ID are equal to the ID of the region that cover its
    """
    # Select n_objects from the objs_map  
    selected_objs_mask, _ = _select_random_objs(objs_map.copy(),
                                                 n_objects,
                                                 remove_obj_border=True)
    
    # Dilate the selected objects to cover all its borders
    selected_objs_mask[selected_objs_mask > 0 ] = 255
    kernel = np.ones((7,7),np.uint8)
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
        forbidden_area[forbidden_area>0] = 255
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
        
        # Apply the copy_move forgery
        if poisson:
            poisson_img, obj_gt_aux, bkg_gt_aux = copy_move_poisson(img,omask, aux_img,background_region)
            if poisson_img is None:
                continue
            aux_img = poisson_img
        else:
            aux_img, obj_gt_aux, bkg_gt_aux = _copy_move(img,omask, aux_img,background_region)
        
        # Update obj_gt_map, bkg_gt_map, and obj_map
        bkg_gt_aux[bkg_gt_aux>0] = gt_index 
        bkg_gt_map= cv2.bitwise_or(bkg_gt_aux,bkg_gt_map)
        obj_gt_aux[obj_gt_aux>0] = gt_index 
        obj_gt_map= cv2.bitwise_or(obj_gt_aux,obj_gt_map)
        objs_map = cv2.bitwise_or(objs_map,obj_gt_map)
        
        # update gt_index
        gt_index +=1
    
    f_img = aux_img
    # Apply post-processing
    f_img = _border_forgery_post_processing(f_img, obj_gt_map)
    
    return f_img, bkg_gt_map, obj_gt_map 



def copy_move_poisson(donor, donor_map, host, host_map,mix=False):
    """
     Perform a copy move forgery copying the donor region (highlighted on donor_map)
    to the host region (highlighted on host_map). 
    
    The copy move will be applied on just one contiguos region
   
    Poisson  use the opencv implementation of Pérez et al. Poisson image editing (2003).
    """
    assert donor.shape[:2] == donor_map.shape[:2] , "Donor image must be the same size of Donor Map"
    assert host.shape[:2] == host_map.shape[:2] , "Host image must be the same size of Host Map"
    
    if len(donor.shape) < 3:
        donor = cv2.cvtColor(donor,cv2.COLOR_GRAY2BGR)
    if len(host.shape) < 3:
        host = cv2.cvtColor(host,cv2.COLOR_GRAY2BGR)
    
        
    # Resize the Donor image to have the same size of the Host
    orig_donor_shape = donor.shape[:2]
    rows, cols = host.shape[:2]
    donor = cv2.resize(donor,(cols,rows))
    donor_map = cv2.resize(donor_map,(cols,rows))
    
    #--------------------------
    # Generate Groundtruth
    #-------------------------
    # Get donor region location
    donor_map[donor_map>0] = 1
    donor_loc = ndi.find_objects(donor_map)[0]
    donor_bb_y, donor_bb_x = _get_start_x_y_from_slice(donor_loc)
    
    # Get Host region location
    host_loc = ndi.find_objects(host_map)[0]
    host_bb_y, host_bb_x = _get_start_x_y_from_slice(host_loc)
    
    # Translate the donor to the position of the host 
    x_shift = host_bb_x - donor_bb_x 
    y_shift = host_bb_y - donor_bb_y 
    
    # Translation Matrix
    M = np.float32([[1, 0, x_shift], [0, 1, y_shift]]) 
    # Get shape of host_map
    rows, cols = host_map.shape[0:2]
    # Donor Translation 
    donor_map_translated = cv2.warpAffine(donor_map.copy(), \
                                  M, dsize=(cols, rows)) 
    
    # Host could be on the edge of the image, so the donor region
    # could appear cropped
    # Assert the shape of the image on the host location
    copy_move_map = cv2.bitwise_and( donor_map_translated, host_map)
    
    # Make the inverse operation to get the real map of  the
    # donor groundtruth region
    M_1 = np.float32([[1, 0, -x_shift], [0, 1, -y_shift]]) 
    donor_gt = cv2.warpAffine(copy_move_map.copy(), M_1, dsize=(cols, rows))
    
    # Get donor region location
    donor_gt[donor_gt>0] = 255
    
    # Get Host region location
    copy_move_map[copy_move_map>0] = 1
    host_point = ndi.center_of_mass(copy_move_map)
    host_point = (int(host_point[1]), int(host_point[0]))
    
    #--------------------------
    # Apply copy_move operation
    #--------------------------
    try:
        if mix:
            f_img = cv2.seamlessClone(donor, host.copy(), donor_gt, host_point, cv2.MIXED_CLONE)
        else:
            f_img = cv2.seamlessClone(donor, host.copy(), donor_gt, host_point, cv2.NORMAL_CLONE)
    except:
        # seamlessClone raised an error because the cloned region 
        # falls outside the limits of the image
        return None,None,None
        
    
    # Resize Donor Gt to original size
    rows, cols = orig_donor_shape
    donor_gt = cv2.resize(donor_gt,(cols,rows))
    
    # Forgery image, Donor's Groundtruth, Host's Groundtruth
    return f_img, donor_gt, np.array(copy_move_map)
    


    

def _translation(im, obj_mask,t_coord,scale=1):
    
    if scale == 0:
        raise ZeroDivisionError("Scale can not be zero.")
    
    
    # Get shape of the image
    rows, cols = im.shape[0:2]

    # Translation matrix
    M = np.float32([[scale, 0, scale*t_coord[0]], [0, scale, scale*t_coord[1]]]) 
    # Execute translation
    result = cv2.warpAffine(im.copy(), M, dsize=(cols, rows)) 
    result_mask = cv2.warpAffine(obj_mask, M, dsize=(cols, rows))
    
    # Avoid overlap between transformated objects
    # If some object is translated on top of  another,
    # we remove the overlap, keeping the object as it was before the
    # transformation.
    # Avoid  overlap on the maks
    result_mask[obj_mask>0] = 0
    # Avoid overlap on the result image
    if len(result.shape) > 2:
        mask_rgb = np.dstack([obj_mask for i in range(result.shape[-1])])
    else:
        mask_rgb = obj_mask
    result = np.where(mask_rgb>0, im, result)
    
    # Since some object could be placed outside the image,
    # We have to refresh the original obj_mask with the inverted
    # transformations
    M_1 = cv2.invertAffineTransform(M)
    refresh_mask = cv2.warpAffine(result_mask, M_1, dsize=(cols, rows))
    
    return result, result_mask, refresh_mask
    
def _rotation(im, obj_mask,rotate,scale):
    
    
    # Get shape of the image
    rows, cols = im.shape[0:2]
    # Translation matrix
    M = cv2.getRotationMatrix2D((cols/2,rows/2),rotate,scale)
    # Execute Rotation
    result_mask = cv2.warpAffine(obj_mask, M, (cols,rows))
    result = cv2.warpAffine(im, M, (cols,rows))
    
    result_mask[obj_mask>0] = 0
    # Avoid overlap on the result image
    if len(result.shape) > 2:
        mask_rgb = np.dstack([obj_mask for i in range(result.shape[-1])])
    else:
        mask_rgb = obj_mask
    result = np.where(mask_rgb>0, im, result)
    
    # Since some object could be placed outside the image,
    # We have to refresh the original obj_mask
    M_1 = cv2.invertAffineTransform(M)#cv2.getRotationMatrix2D((cols/2,rows/2),-rotate,1/scale)
    refresh_mask = cv2.warpAffine(result_mask, M_1, dsize=(cols,rows))
    
    return result, result_mask, refresh_mask

def _mirrow(im, obj_mask):
    
    # Flip image and mask
    result_mask = cv2.flip(obj_mask, 1)
    result = cv2.flip(im, 1)
    
    # Remove overlaped objects
    result_mask[obj_mask>0] = 0
    # Avoid overlap on the result image
    if len(result.shape) > 2:
        mask_rgb = np.dstack([obj_mask for i in range(result.shape[-1])])
    else:
        mask_rgb = obj_mask
    result = np.where(mask_rgb>0, im, result)
    
    # Since some object could be placed on top of each other
    # We have to refresh the mask
    refresh_mask = cv2.flip(result_mask, 1)
    
    return result, result_mask, refresh_mask

def _set_transformations(im, obj_mask, t_types):
    """
    Since the object after the operation could
    have some of its content falling out the image,
    we do not allow the operation to be done at once.
    Otherwise, the result image could have visible distortion,
    that an human easily could point.
    """
    
    forgery_mask = np.zeros_like(obj_mask)
    
    # Perform flip on image
    if t_types['flip']:
            im, forgery_mask, obj_mask = _mirrow(im, obj_mask)
            
     # Perform translation
    elif isinstance(t_types['translation'], tuple):
        # If scale is included on transformation
        if isinstance(t_types['scale'], int)\
            or isinstance(t_types['scale'], float):
            scale = t_types['scale']
        else:
            scale = 1
            
        im, forgery_mask, obj_mask = _translation(im, obj_mask,t_types['translation'],scale)
        
    # Perform Rotation
    elif isinstance(t_types['rotation'], int):
        # If scale is included on transformation
        # do scaling along with rotation
        if isinstance(t_types['scale'], int)\
            or isinstance(t_types['scale'], float):
            scale = t_types['scale']
        else:
            scale = 1
        im, forgery_mask, obj_mask = _rotation(im, obj_mask, t_types['rotation'], scale)
    
    else:
        raise NotImplementedError
        
    
   
    return im, forgery_mask, obj_mask
    

def copy_move_forgery(img, objs_map, t_types, n_objects=1):
    """
    Perform a copy move forgery on 'img' following the transformation pass 
    on the input dictionary t_types.
    
    Parameters
    ---------
    img: opencv python object 
    objs_map:opencv python object 
        binary image containing the segmentated version of the foreground
    n_objects: int
        Number of objects to perform the copy_move
    t_types: dict - Transformation Types
        {
        'translation': Tuple(<int>,<int>) in [im.shape[0],im.shape[1]],
        'rotation': <int>  in [0,360],
        'scale' : <float>,
        'flip': <bool>,
        }
        Note that 'flips', 'translation' and 'rotation' are individually
        applied. To avoid many overlaps between the objects, we not allow those 
        transformatios to be applied together on the same image.
        
    Result
    ------
    f_img: cv2 image obj
        Image with copy move forgery
    f_mask: cv2 image obj
        Grayscale image forgery, where each grayscale represent where the forgery object id appears,
        yet 0 that represents the background
    obj_mask: cv2 image obj
        Grayscale image forgery, where each grayscale represent where the origianl object id 
        that were selected during the forgery process.
        0 represents the background
        
    """
    
    #Save orig image
    orig_im = img.copy()
    
    # Create defaultdict with the input dict
    t_types = defaultdict(lambda: None,t_types)

    # select random objects form 
    obj_mask, _ = _select_random_objs(objs_map, n_objects,remove_obj_border=True)
    
    obj_mask[obj_mask >0] = 1
    obj_mask,_  = ndi.label(obj_mask)
    obj_mask = obj_mask.astype(np.uint8)
    #------------------------#
    # Perform Transformations
    #------------------------#
   
    forgery_img, f_mask, obj_mask = _set_transformations(img, obj_mask, t_types)
   
        
    # Get a PIL version of the original img
    img_pil = Image.fromarray(orig_im)
    # Mask and transformated img
    mask = f_mask.copy()
    mask[mask>0] = 255
    mask = Image.fromarray(mask)
    
    forgery_pil = Image.fromarray(forgery_img)

    # Make Forgery
    img_pil.paste(forgery_pil, (0, 0),mask)
    
    f_img = np.array(img_pil)
    
    # Apply post-processing
    f_img = _border_forgery_post_processing(f_img, f_mask)
    
    return f_img, f_mask, obj_mask


def _border_forgery_post_processing(f_img, f_mask):
    """
    Smooth the edges of each object to avoid a hard glitch
    between the forgery objec and the image background.
    """
    
    # Get the border of each object
    kernel = np.ones((3,3),np.uint8)
    erode_mask = cv2.erode(f_mask,kernel,iterations = 1)
    dilate_mask = cv2.dilate(f_mask,kernel,iterations = 1)
    xor_mask = cv2.bitwise_xor(erode_mask,dilate_mask)
    
    # Find the Bounding box of each element
    xor_mask[xor_mask>0]=1
    xor_mask,_ = ndi.label(xor_mask)
    bboxs = ndi.find_objects(xor_mask.astype(np.int32)) 
    
    # Apply a GaussianBlur on its borders
    f_blur = f_img.copy()
    for bbox in bboxs:
        f_blur[bbox] = cv2.GaussianBlur(f_blur[bbox],(3,3),0)
        
    # Apply the blur border on the f_img    
    xor_mask[xor_mask>0] = 255
    xor_mask = xor_mask.astype(np.uint8)
    xor_mask = cv2.bitwise_not(xor_mask).astype(bool)
    if len(xor_mask.shape) != len(f_blur.shape):
        xor_mask = np.dstack([xor_mask for i in range(f_blur.shape[-1])])
    np.putmask(f_blur, xor_mask, f_img)
    
    return f_blur