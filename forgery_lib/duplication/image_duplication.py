"""
Image Duplication Forgery
This is a set of functions that implements types of duplication forgeries
on digital images for reuse manipulation.

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
from .overlap import _perform_crop

def simple_copy(img, operations):
    """
    Perfom a simple copy of the entery input image addressing any operation from the dict 
    operations.
    
    Parameters
    ---------
    img: cv2 source img
    operations:dict
        Operations applied on mask1 contains:
        {
            flip: <boolean>
            rotation: <False|| None || CV2 rotation flag>
                -> [ 'ROTATE_180' , 'ROTATE_90_CLOCKWISE', 
                            'ROTATE_90_COUNTERCLOCKWISE']
            brightness: <int>
            }
    Return
    ------
    src_img: cv2 img obj
        image region of 'im' using operations from 'operation1'
    src_mask1: Binary image 
        represent the regions from the source image that was manipulated
        before the manipulation
    
    forgery_img: cv2 img obj
        Forgery version of the src_img using 'operations'
    forgery_gt: Binary image
        represent the regions from the source image that was manipulated
        after the manipulation
        
    Example
    -------
    op = {
            'flip': False,
            'rotation': ROTATE_180,
            'brightness': 20
            }

    src_img, src_gt , forgery_img, forgery_gt = simple_copy(img, op)
    """
    # Create defaultdict with the input dict
    operations = defaultdict(lambda: None,operations)
    operations['crop_type'] = None

    # operation done on mask
    copy, mask = _perform_crop(img,operations)
    
    return img,mask.copy(),\
        copy,mask.copy()



 

