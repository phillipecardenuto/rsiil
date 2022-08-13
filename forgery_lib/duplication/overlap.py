"""
Overlap Forgeries
This is a set of functions that implements types of overlap forgeries
on digital images.

-------------------------------------
Author: Phillipe Cardenuto
email: phillipe.cardenuto@gmail.com
--------------------------------------
Version:                           1.0
Version Date:               July, 2020

This is file is part of Forgery Library
"""

import cv2
from scipy import ndimage as ndi
import numpy as np
from collections import defaultdict
from ..utils import increase_brightness


def _undo_rot_flip(mask, rotation,flip):
       
# If Flip was applied, correct gt
    if flip:
        mask = cv2.flip(mask,1)

    # If Rotation was applied, correct gt
    elif rotation in ["ROTATE_180","ROTATE_90_CLOCKWISE","ROTATE_90_COUNTERCLOCKWISE"]\
               and not isinstance(rotation,bool):
        if rotation == "ROTATE_180":
            undo_rot = cv2.ROTATE_180
        elif rotation == "ROTATE_90_CLOCKWISE":
            undo_rot = cv2.ROTATE_90_COUNTERCLOCKWISE
        elif rotation == "ROTATE_90_COUNTERCLOCKWISE":
            undo_rot = cv2.ROTATE_90_CLOCKWISE
        # Mask Correction
        mask = cv2.rotate(mask,undo_rot)
            
    return mask

def _do_rot_flip(im, rotation,flip):
    # Flip
    if flip:
        im = cv2.flip(im,1)
        
    # Rotation
    elif rotation in ["ROTATE_180","ROTATE_90_CLOCKWISE","ROTATE_90_COUNTERCLOCKWISE"]\
        and not isinstance(rotation,bool):
        cv2_flag = eval(f"cv2.{rotation}")
        im = cv2.rotate(im,cv2_flag)
    
    return im
    

#Perform crop with ground truth on 
def _perform_crop(im, operations):
    """
    Perform crop and  generate a binary map of the region croped
    related to the original image 'im'
    
    Parameters
    ---------
    im: cv2 img
    
    operatios:dict
      Operation applied on im before croping the image 'im'
            {
            flip:boolean,
            rotation: None or CV2 rotation flag
                 -> [ 'cv2.ROTATE_180' , 'cv2.ROTATE_90_CLOCKWISE', 
                            'cv2.ROTATE_90_COUNTERCLOCKWISE'],
            crop_type:str in 
                    -> ['top_right', 'top_left', 'center',None]
            }
    """
    #TODO: Reject false operations
    
    # Get width and height
    width, height = im.shape[0], im.shape[1]
    
    # Perform rotation and flip
    im = _do_rot_flip(im,operations['rotation'],operations['flip'])
    # Ground Truth img
    gt = np.zeros_like(im)
    
    # Perform crop operation
    if operations['crop_type'] == 'top_left':
        crop  = im[0:width//2, 0:height//2] 
        gt[0:width//2, 0:height//2] = 255
    elif operations['crop_type'] == 'top_right':
        crop  = im[0:(3*width)//5,(2*height)//5: ] 
        gt[0:(3*width)//5,(2*height)//5:]  = 255
    elif operations['crop_type'] == 'center':
        crop  = im[(1*width)//5:(7*width)//10, (1*height)//5:(7*height)//10]
        gt[(1*width)//5:(7*width)//10, (1*height)//5:(7*height)//10] = 255
    else:
        crop  = im
        gt = gt + 255
        
    # Change Brightness if requested
    if operations['brightness']:
        crop = increase_brightness(crop, operations['brightness'])
    
    # Fix groundtruth
    gt = _undo_rot_flip(gt, operations['rotation'], operations['flip'])
            
    return crop, gt

def _crop_mask(mask, crop_type, width, height):
    """
    Width and height are relative to the original image and not to the crop
    """
    
    # Get width and height
    #width,height = mask.shape[:2]
    
    if crop_type == 'top_left':
        return mask[0:width//2, 0:height//2]
    
    elif crop_type == 'top_right':
        return mask[0:(3*height)//5, (2*width)//5: ]
    
    elif crop_type == 'center':
        return mask[(1*width)//5:(7*width)//10, (1*height)//5:(7*height)//10]
    
    else:
        return mask
        

def overlap_forgery(img,operation1,operation2):
    """
    From a source image 'img' perform two types of operations aiming to create 2 cropped versions
    from 'img'.
    Flip and rotaiton are not allowed at the same time. If both were declared just the flip will be done.
    
    Parameters
    ---------
    img: cv2 img
    operation1:dict
        Operations applied on mask1 contains:
            {
            flip: <boolean>
            rotation: <False|| None || CV2 rotation flag>
                 -> [ 'ROTATE_180' , 'ROTATE_90_CLOCKWISE', 
                            'ROTATE_90_COUNTERCLOCKWISE']
            crop_type:<str>  
                    -> ['top_right', 'top_left', 'center']
                    
            brightness: <int>
            }
    operation2: dict
        Operations applied on mask2 contains:
        {
            flip: <boolean>
            rotation: <False|| None || CV2 rotation flag>
                 -> [ 'ROTATE_180' , 'ROTATE_90_CLOCKWISE', 
                            'ROTATE_90_COUNTERCLOCKWISE']
            crop_type:<str>  
                    -> ['top_right', 'top_left', 'center']
                    
            brightness: <int>
            }
    Return
    ------
    crop1: cv2 img obj
        Cropped image region of 'img' using operations from 'operation1'
    overlap_mask1: Binary image 
        represent the portion of the overlap of the mask1 to mask2
    
    crop2: cv2 img obj
        Cropped image region of 'img' using operations from 'operation2'
    overlap_mask2: Binary image
        represent the portion of the overlap of the mask2 to mask1
        
    Example
    -------
    op1 = {
            'crop_type': 'top_left',
            'flip': False,
            'rotation': False}
    op2 = {
            'crop_type': 'top_right',
            'flip': False,
            'rotation': cv2.ROTATE_180 
            }
    crop1, gt_1 ,crop2, gt_2 = overlap_forgery(img, op1, op2)
    """
    
    # Create defaultdict with the input dict
    operation1 = defaultdict(lambda: None,operation1)
    operation2 = defaultdict(lambda: None,operation2)

    # redo operation done on mask
    crop1, mask1 = _perform_crop(img,operation1)
    crop2, mask2 = _perform_crop(img,operation2)
    
    # Apply an 'and' bitwise  operation to locate the portion of the image that overlaps
    resultant_overlap = cv2.bitwise_and(mask1,mask2)
    
    mask1 = _do_rot_flip(resultant_overlap,operation1['rotation'],operation1['flip'])
    mask2 = _do_rot_flip(resultant_overlap,operation2['rotation'],operation2['flip'])
    
    # Get the croped_version of the masks
    croped_mask1 = _crop_mask(mask1,operation1['crop_type'], img.shape[0], img.shape[1])
    croped_mask2 = _crop_mask(mask2,operation2['crop_type'], img.shape[0], img.shape[1])
    
    return crop1,croped_mask1,\
           crop2,croped_mask2
