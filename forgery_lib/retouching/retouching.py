"""
Retouching Forgeries
This is a set of functions that implements 
retouching forgeries.

-------------------------------------
Author: Phillipe Cardenuto
email: phillipe.cardenuto@gmail.com
--------------------------------------
Version:                           1.1
Version Date:               Aug , 2020

This is file is part of Forgery Library
"""

import cv2
from ..utils import  _select_random_objs
from scipy import ndimage as ndi
import numpy as np
from PIL import Image

def retouching_blurring(img, objs_map, n_objects=1, 
                        rand=True, blur_level=4,
                        gaussian=False):
    """
    Apply hard retouching blurring method on image 
    aiming to obfuscate some objects
    
    Parameters
    ---------
    img: <np.array> Image input
    objs_map: <np.array> Map each pixel for an Object ID,
            or highlight the object (255) from the background (0)
    n_objects: <int> 
        number of objects that the retouching will be applied
    rand: optional <bool> = True
        Select object at random
    blur_level: optional <int> = 4
        How blurry the result region would be transformed
    gaussian: optional <bool> = False
        If True apply a gausian filter
        
    Return
    ------
    f_img: opencv obj
        Image with retouching forgery
    gt: opencv obj
        Binary image highlighting the regions that retouching was applied
        255 -> Forgery regions
        0   -> normal
    """
    
    
    # Select Random objects 
    if rand:
        objs_map, obj_location = _select_random_objs(objs_map, n_objects,remove_obj_border=False)
        
    # Get the location of all objs on image
    else:
        objs_map[objs_map>0] = 255
        obj_labels, total_obj = ndi.label(objs_map)
        obj_location = ndi.find_objects(obj_labels)
        
    
    # Apply blurring function
    # Ref: https://docs.opencv.org/master/d4/dbd/tutorial_filter_2d.html
    retouching_img = img.copy()
    retouching_map = np.zeros_like(objs_map)
    
    # Define kernel
    kernel_size = 11
    kernel = np.ones((kernel_size, kernel_size), dtype=np.float32)
    kernel /= (kernel_size * kernel_size)
    
    for index, obj in enumerate(obj_location):
        # Hard Blurring
        for i in range(blur_level):
            if gaussian:
                retouching_img[obj] =  cv2.GaussianBlur(retouching_img[obj], (11,11),0)
            else: 
                retouching_img[obj] = cv2.filter2D(retouching_img[obj], -1, kernel)
        
        # Insert modified region on retouching map
        retouching_map[obj] = 255
        
    retouching_img = retouching_img.astype(np.uint8)
    # Apply the forgery only on the exact shape of the object
    # generating a ground truth for the forgeries
    return _apply_retouching_forgery(img, objs_map,
                          retouching_img, retouching_map)    

    
def retouching_brightness_contrast(img, objs_map, n_objects=1, alpha=None, beta=None):
    """
    Apply hard retouching brightness and contrast method on image, 
    aiming to obfuscate some objects
    Alpha and Beta represent the contrast and brightness respectively
    Alpha and Beta are list that contains n_objs. Each respective value of
        those lists will be applied to an object.
        
    Parameters
    ---------
    img: <np.array> Image input
    objs_map: <np.array> Map each pixel for an Object ID,
            or highlight the object (255) from the background (0)
    n_objects: <int> 
        number of objects that the retouching will be applied
        
   alpha: list of floats [<float>] [0,5]
       Alpha value for the contrast
       If alpha is None, random value will be generated
   beta: list of ints [<int>] [0,255]
       Beta value for brightness
       If beta is None, random value will be generated
    
    Return
    ------
    f_img: opencv obj
        Image with retouching forgery
    gt: opencv obj
        Binary image highlighting the regions that retouching was applied
        255 -> Forgery regions
        0   -> normal
    """
    
    # If alpha was not defined, get a random value 
    if alpha is None:
        #numpy.random.rand() * (max - min) + min --> Random value between [max,min]
        # [0,2]
        alpha = [np.random.random()*(2) for i in range(n_objects)]
    assert type(alpha) is list, "Alpha must be a list"
    
    while len(alpha) < n_objects:
        alpha.append(np.random.random()*(2))
    
    
    # If beta was not defined, get a random value 
    if beta is None:
        beta = [ np.random.randint(20,100) for i in range(n_objects)]
    assert type(beta) is list, "Beta must be a list"
    
    while len(beta) < n_objects:
        beta.append(np.random.randint(20,100))
    
    # Select Random objects 
    objs_map, obj_location = _select_random_objs(objs_map, n_objects,remove_obj_border=False)
    
    # Apply contrast function
    #Ref https://docs.opencv.org/3.4/d3/dc1/tutorial_basic_linear_transform.html
    retouching_img = img.copy().astype(np.uint32)
    retouching_map = np.zeros_like(objs_map)
    
    for index, obj in enumerate(obj_location):
        # Modify Contrast 
        retouching_img[obj] = (alpha[index] * retouching_img[obj] ) + beta[index]
        # Assert that the img fit in np.uint8
        retouching_img[retouching_img > 255] = 255
        retouching_img[retouching_img < 0] = 0
        
        # Insert modified region on retouching map
        retouching_map[obj] = 255
        
    retouching_img = retouching_img.astype(np.uint8)
    # Apply the forgery only on the exact shape of the object
    # generating a ground truth for the forgeries
    return _apply_retouching_forgery(img, objs_map,
                          retouching_img, retouching_map)    

def retouching_brightness(img, objs_map, n_objects=1, bright_level=None):
    """
    Apply hard retouching brightness method on object from an image.
    
    Parameters
    ---------
    img: <np.array> Image input
    objs_map: <np.array> Map each pixel for an Object ID,
            or highlight the object (255) from the background (0)
    n_objects: <int> 
        number of objects that the retouching will be applied
   bright_level: (optional) list of int [<int>] [-255,255] 
       Level of brightness of each obj. If the size of the list is smaller than the n_objexts, a random value will
       be generated for each missing object's brightness.
    
    Return
    ------
    f_img: opencv obj
        Image with retouching forgery
    gt: opencv obj
        Binary image highlighting the regions that retouching was applied
        255 -> Forgery regions
        0   -> normal
    """
    
    # If bright_level was not defined, get a random values
    if bright_level is None:
        bright_level = [ np.random.randint(-255,255) for i in range(n_objects)]
    else:
        assert type(bright_level) is list, "bright_level must be a list"
        bright_level = np.array(bright_level).astype("int8")


    
    while len(bright_level) < n_objects:
        bright_level.append(np.random.randint(-255,255))
        
    # Select Random objects 
    objs_map, obj_location = _select_random_objs(objs_map, n_objects, remove_obj_border=False)
    
    # Apply Brightness Adjust
    retouching_img = img.copy().astype(np.int32)
    retouching_map = np.zeros_like(objs_map)
    
    for index, obj in enumerate(obj_location):
        # Modify Contrast 
        
        retouching_img[obj]+= bright_level[index] 
        # Assert that the img fit in np.uint8
        retouching_img[retouching_img > 255] = 255
        retouching_img[retouching_img < 0] = 0
        
        # Insert modified region on retouching map
        retouching_map[obj] = 255
    
    retouching_img = retouching_img.astype(np.uint8)
    
    # Apply the forgery only on the exact shape of the object
    # generating a ground truth for the forgeries
    return _apply_retouching_forgery(img, objs_map,
                          retouching_img, retouching_map)    
     

def retouching_contrast(img, objs_map, n_objects=1, alpha=None):
    """
    Apply hard retouching contrast method on image, 
    aiming to obfuscate some objects
    
    Parameters
    ---------
    img: <np.array> Image input
    objs_map: <np.array> Map each pixel for an Object ID,
            or highlight the object (255) from the background (0)
    n_objects: <int> 
        number of objects that the retouching will be applied
   alpha: <float> [0,3] 
       Alpha value for the contrast
    
    Return
    ------
    f_img: opencv obj
        Image with retouching forgery
    gt: opencv obj
        Binary image highlighting the regions that retouching was applied
        255 -> Forgery regions
        0   -> normal
    """
    
    # If alpha was not defined, get a random value 
    if alpha is None:
        # [0,3]
        alpha = np.random.random()*(3)
        
    # Select Random objects 
    objs_map, obj_location = _select_random_objs(objs_map, n_objects, remove_obj_border=False)
    
    # Apply contrast function
    #Ref https://huningxin.github.io/opencv_docs/d3/dc1/tutorial_basic_linear_transform.html
    retouching_img = img.copy().astype(np.uint32)
    retouching_map = np.zeros_like(objs_map)
    
    for obj in obj_location:
        # Modify Contrast 
        retouching_img[obj] = alpha * retouching_img[obj] 
        # Assert that the img fit in np.uint8
        retouching_img[retouching_img > 255] = 255
        retouching_img[retouching_img < 0] = 0
        
        # Insert modified region on retouching map
        retouching_map[obj] = 255
        
    retouching_img = retouching_img.astype(np.uint8)
    # Apply the forgery only on the exact shape of the object
    # generating a ground truth for the forgeries
    return _apply_retouching_forgery(img, objs_map,
                          retouching_img, retouching_map)    
     

def _apply_retouching_forgery(orig_img, obj_map, forgery_img, forgery_map):

    # Dilate the obj_map to cover the object border
    kernel = np.ones((7,7),np.uint8)
    obj_map = cv2.dilate(obj_map,kernel,iterations = 1)
    obj_map[obj_map>0] = 255
    
    # Applied the forgery only on selected regions from 
    # the forgery map ( Foregey map should be rectangular)
    # The next line give the shape of the obj
    # to the foregery map
    forgery_map = cv2.bitwise_and(forgery_map, obj_map)
    
    # Mask and transformated img
    forgery_map[forgery_map>0] = 255
    # Use of PIL to create the forgery
    forgery_map_pil= Image.fromarray(forgery_map)
    
    
    # Open the final forgery image as PIL image
    result_img_pil = Image.fromarray(orig_img)
    # Open the foregery image as PIL image
    forgery_img_pil = Image.fromarray(forgery_img)
    # Make Forgery
    result_img_pil.paste(forgery_img_pil, (0, 0),forgery_map_pil)
    
    # Pil to numpy
    result_foregery_img = np.array(result_img_pil)
    result_gt = np.array(forgery_map_pil)
    
    return result_foregery_img, result_gt