"""
Utils
This is a set of functions that is usefull for
debugging or auxiliar for the forgery-lib

-------------------------------------
Author: Phillipe Cardenuto
email: phillipe.cardenuto@gmail.com
--------------------------------------
Version:                           1.0
Version Date:               July, 2020

This is file is part of Forgery Library
"""
import cv2
from PIL import Image
import numpy as np
from scipy import ndimage as ndi

def increase_brightness(img, value):
    """
    Reference: https://stackoverflow.com/a/47427398
    """
    GRAY=False
    ALPHA=False
    if len(img.shape) == 2:
        GRAY=True
        img = cv2.cvtColor(img, cv2.COLOR_GRAY2RGB)
    elif img.shape[-1] > 3:
        ALPHA = True
        alpha = img[:,:,-1]
        img = img[:,:,:3]
    hsv = cv2.cvtColor(img, cv2.COLOR_RGB2HSV)
    h, s, v = cv2.split(hsv)

    lim = 255 - value
    v[v > lim] = 255
    v[v <= lim] += value

    final_hsv = cv2.merge((h, s, v))
    img = cv2.cvtColor(final_hsv, cv2.COLOR_HSV2RGB)
    
    if GRAY:
        img = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
    if ALPHA:
        r,g,b = cv2.split(img)
        img = cv2.merge((r,g,b,alpha))
    return img


def exposure(im, alpha):
    """
    Apply a exposure on im.
    This method could be applied on image to visually
    check if any forgery was applied on a image.

    Parameters
    ----------
    im: opencv object
        Image that the exposure will be applied
    alpha:
        Level of exposure to be applied

    Return
    ------
    e_img: img np.float 32
    """

    # Exposure image before method
    e_img = im.copy().astype(np.float32)

    # Image that allow values higher then 255
    overflow = e_img.copy().astype(np.float32)
    overflow = overflow *(2**alpha)

    # e_img could not allow overflow values higher then 255
    e_img = np.where((e_img*(2**alpha) > 255), np.ones_like(e_img)*255, overflow)
    return e_img


# Debug panel
def display_panel(imgs, display_type='row'):
    """
    Easy way to create a panel to check any image from a imgs list

    Parameters
    ----------
    imgs: list
        List of image on opencv format (numpy)

    display_type: str <'row'| 'col'>
        Type of the order to display the images from the 'imgs' list.
        'row': Row based sort
                1 2 3
                4 5 6
        'col': Column based sort
                1 3 5
                2 4 6
    Return
    ------
    panel: PIL object
        A PIL Image object with the imgs inputed in a panel format
    """
    # Total imgs
    total_imgs = len(imgs)

    # Get min shape
    min_shape = sorted( [(np.sum(i.shape[:2]), i.shape[:2] ) for i in imgs])[0][1]
    rows,cols = min_shape
    # Adding border to imgsa
    bordersize = 2

    imgs_aux = [ cv2.copyMakeBorder(
                            cv2.resize(im,(cols,rows)), # Resize img
                            top=bordersize,
                            bottom=bordersize,
                            left=bordersize,
                            right=bordersize,
                            borderType=cv2.BORDER_CONSTANT,
                            value=[69, 87, 96]

                    ) for im in imgs]
    # number of img in vertical
    nv_imgs = int(np.floor(np.sqrt(total_imgs)))
    # number of img in horizontal
    nh_imgs = int(np.ceil(total_imgs/nv_imgs))

    # Empty imgs
    empty =  nv_imgs*nh_imgs - total_imgs


    # Image Panel on array format
    img_ar = imgs_aux[0]

    # Insert missing imgs to complete a rectangle imgs
    if empty > 0:
        for i in range(empty):
            ii = np.zeros_like(img_ar)
            ii[:] = [69, 87, 96]
            imgs_aux += [ii]

    # Display  panel using row based sort
    if display_type == 'row':
        # Concatenate all imgs
        for i in range(1,len(imgs_aux)):
            img_ar = np.concatenate([img_ar,imgs_aux[i]],axis=1)

        # Reshape properly to (nh_imgs nv_imgs)
        r = img_ar.reshape((img_ar.shape[0],nv_imgs,-1)+img_ar.shape[2:])
        if len(r.shape) > 3:
            panel = np.vstack([r[:,i,:,:] for i in range(nv_imgs)])
        else:
            panel = np.vstack([r[:,i,:] for i in range(nv_imgs)])



    # Display  panel using column based sort
    else:
        # Concatenate all imgs
        for i in range(1,len(imgs_aux)):
            img_ar = np.concatenate([img_ar,imgs_aux[i]],axis=0)

        # Reshape properly to (nh_imgs nv_imgs)
        r = img_ar.reshape((nh_imgs,-1)+img_ar.shape[1:])
        panel =  np.hstack(r)

    return  panel


def _remove_border_objects(obj_map):
    """
    Remove objects that touch the border
    Those objects probably are cut, and will easy
    distinguish from the rest
    """

    # Label each object from object map
    obj_labels, n_obj = ndi.label(obj_map)

    discard = list([
                    np.unique(obj_labels[0,:]), # Top Border
                    np.unique(obj_labels[:,-1]), # Right Border
                    np.unique(obj_labels[-1,:]),# Bottom Border
                    np.unique(obj_labels[:,0]) # Left Border
                    ])

    discard = np.unique(np.hstack(discard))

    for i in discard:

        if i != 0:
            obj_labels[obj_labels == i] = 0
    # Reordering ids
    obj_labels[obj_labels > 0] = 255

    obj_labels, n_obj = ndi.label(obj_labels)

    return obj_labels, n_obj

def _get_start_x_y_from_slice(sl):
    y = sl[0].start
    x = sl[1].start
    return y,x

def _find_background(gt,obj_bb):
    """
    Find a region in gt that fits the obj_bb
    Parameters
    ----------
    gt: map with the objects (foreground) with non-zero values
    obj_bb: object bounding box
    """
    bb_h, bb_w = obj_bb
    gt_h, gt_w = gt.shape[0:2]
    regions_available = []

    for y in range(0,gt_h,bb_h):
        for x in range(0, gt_w, bb_w):

            if (y + bb_h) <= gt_h  and \
               (x + bb_w) <= gt_w:

                tile = gt[y:y+bb_h,x:x+bb_w]

                if not tile.any():
                    reg = (slice(y,y+bb_h), slice(x,x+bb_w))
                    regions_available.append(reg)

    if regions_available:
        gt_tile = np.zeros_like(gt)
        random_region = np.random.choice(len(regions_available))
        region = regions_available[random_region]
        gt_tile[region] = 1
        return gt_tile

    # None region found
    return None




def _select_random_objs(obj_map, n_objs=1, remove_obj_border=True):
    """
    Select at random  the bbox of n_objs from the image based
    on the obj_map.

    Parameter
    ---------
    obj_map:<np.array>
        Numpy array mapping each location of each
        object from img
    n_objs: <int>
        Number of objects to select from the obj_map
    remove_obj_border: <bool>
        If true remove all objects that touchs the border of the obj_map
    Return
    ------
    The slice of min(n_obj,total_obj) selected in the obj_map

    """
    # Preprocessing Forgery
    # Remove objects on image border and
    # give each object an id
    if remove_obj_border:
        obj_labels, total_obj = _remove_border_objects(obj_map)
    else:
        obj_map[obj_map>0] = 255
        obj_labels, total_obj = ndi.label(obj_map)

    assert total_obj > 0 , f"Number of object in map is 0"

    # Randomly selects some objects
    obj_id = []
    ids = list(np.arange(1,total_obj))
    while len(obj_id) < total_obj and len(obj_id) < n_objs and len(ids) > 0:
        selected_id = np.random.choice(ids)
        obj_id.append(selected_id)
        ids.remove(selected_id)
    # Make a mask image containing only the
    # selected object
    obj_mask = np.isin(obj_labels, obj_id)

    # Since the obj_id could be more then 255
    # Assert that obj_mask is uint8
    obj_mask[obj_mask>0] = 255
    obj_mask = obj_mask.astype(np.uint8)

    # Dilate the regions of each object, to make the applied operation
    # more accurate since the object map could be smaller than the object.
    kernel = np.ones((11,11),np.uint8)
    obj_dilate = cv2.dilate(obj_mask,kernel,iterations = 4)

    # Get the bbox of each selected object (slice form)
    obj, _ = ndi.label(obj_dilate)
    objs_bbox =  ndi.find_objects(obj)


    return obj_mask,objs_bbox

