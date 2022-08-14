import cv2
from scipy import ndimage
from PIL import Image
import os
import pandas as pd
import numpy as np
from tqdm import tqdm

###################
#     METRICS    #
##################

def F1(detected, groundtruth):
    """
    Implementation of f1-score
    considering a detection map and a ground-truth binary mask
    """
    
    activation_map = detected.copy()
    gt = groundtruth.copy()
    activation_map[activation_map>0] = 1
    gt[gt>0] = 1
    # True positives
    tp = cv2.bitwise_and(activation_map, gt).sum()
    
    # False positives + False Negatives
    fp_fn = cv2.bitwise_xor(activation_map,gt).sum()
    
    if tp == 0:
        return 0
    
    f1_score = 2*tp / (2*tp + fp_fn)
    
    return f1_score


def precision(detected, groundtruth):
    """
    Implementation of precision
    considering a detection map and a ground-truth binary mask
    """
    
    activation_map = detected.copy()
    gt = groundtruth.copy()
    activation_map[activation_map>0] = 1
    gt[gt>0] = 1
    # True positives
    tp = cv2.bitwise_and(activation_map, gt).sum()
    
    # False Negatives
    not_gt = cv2.bitwise_not(gt)
    not_activation_map = cv2.bitwise_not(activation_map)
    fn = cv2.bitwise_and(not_activation_map, not_gt).sum()
    
     # Calculate False Positives from XOR map
    # XOR = False positives + False Negatives
    XOR = cv2.bitwise_xor(activation_map,gt)
    fp =  cv2.bitwise_and(activation_map,XOR).sum()
    
    if tp == 0:
        return 0
    
    p = tp / (tp + fp)
    return p
    
def F1_CTP(detected, groundtruth, check_cc=True):
    
    """
    This is the implementation of F1 using consistent truth positive (CTP).
    If you set check_cc=True, you are using F1 with CTP, instead of TP; otherwise,
    the regular true positive will be used.
    Using CTP you are checking if both (object and its copy) were detected. If only one 
    of them are detected, the metric will exclude its region from the detected map.
    You might also input a detected map with different object labeling them with unique ID for
    each forgery (object and copy)
    
    Both inputs must have the objects already labeled by an ID in [1, inf].
    The IDs from the detection map and the gt do not need to match.
    
    Parameters
    ---------
    detected: numpy array with the ids of the objects activated as suspicious
    groundtruth:  Ground truth map with objects involved in the forgery labeled with an ID
    check_cc: Verify whether the connected components from the activation map match with any object from gt
    Example:
            detected = [
                        0 2 2 0 1 
                        0 2 2 0 1
                        0 0 0 0 0
                        1 0 0 2 2
                        1 0 0 2 2
                       ]
                               
            groundtruth = [
                        0 1 1 0 0 
                        0 1 1 0 0
                        0 0 0 0 0
                        0 0 0 1 1
                        0 0 0 1 1
                        ]
    the object labeled as 1 on gt was copied from the bottom to the top of the image, so those object will care the same IDS on the ground truth, but they don't need to be the same
    in the detected regions.
    """
    
    
    # Groundtruth is not pointing to any issue in the image
    if groundtruth.any() == False or check_cc == False:
        return F1(detected, groundtruth)
    
    # Make a copy of the imgs to prevent any modification on the real data
    detection_map = detected.copy()
    gt = groundtruth.copy()
        
    detected_ids = np.unique(detection_map)
    detected_ids = np.delete(detected_ids, np.where(detected_ids==0)) # Remove background id
    
    if len(detected_ids) == 0:
        return 0
    
    # Create consistent detected map
    consistent_detected_map = get_consistent_regions(detection_map, gt)    
    consistent_detected_map[consistent_detected_map>0] = 1
    # Create inconsistent detected map by removing the consistent region from the detection_map
    inconsistent_detected_map = detection_map.copy()
    inconsistent_detected_map[consistent_detected_map>0] = 0
    inconsistent_detected_map[inconsistent_detected_map>0] = 1

    # Calculates F1 score
    TP = np.sum(consistent_detected_map)
    FP = np.sum(inconsistent_detected_map)
    gt[gt>0] = 1    
    consistent_detected_map[consistent_detected_map>0] = 1
    FN = np.sum(gt - consistent_detected_map)

    
    
    F1_score = 2*TP / (2*TP + FP + FN)
    
    return F1_score


def precision_CTP(detected, groundtruth, check_cc=True):
    """
    Precision using CTP.  If you set check_cc=True, you are using Precision with CTP, instead of TP; otherwise,
    the regular true positive will be used.
    Using CTP you are checking if both (object and its copy) were detected. If just one 
    of them are detected, the metric will exclude its region from the detected map.
    You might also input a detected map with different object labeling them with unique ID for
    each forgery (object and copy).
    
    Both inputs must have the objects already labeled by an ID in [1, inf].
    The IDs from the detection map and the gt do not need to match.
    
    Parameters
    ---------
    detected: numpy array with the ids of the objects activated as suspicious
    groundtruth:  Ground truth map with objects involved in the forgery labeled with an ID
    check_cc: Verify whether the connected components from the activation map match with any object from gt
    Example:
            detected = [
                        0 2 2 0 1 
                        0 2 2 0 1
                        0 0 0 0 0
                        1 0 0 2 2
                        1 0 0 2 2
                       ]
                               
            groundtruth = [
                            0 1 1 0 0 
                            0 1 1 0 0
                            0 0 0 0 0
                            0 0 0 1 1
                            0 0 0 1 1
                          ]
    the object labeled as 1 on gt was copied from the bottom to the top of the image, so those object will care the same IDS on the ground truth, but they don't need to be the same
    in the activated regions.
    
    """
    
    # Groundtruth is not pointing to any issue in the image
    if groundtruth.any() == False or check_cc == False:
        return F1(detected, groundtruth)
    
    # Make a copy of the imgs to prevent any modification on the real data
    detection_map = detected.copy()
    gt = groundtruth.copy()
        
    detected_ids = np.unique(detection_map)
    detected_ids = np.delete(detected_ids, np.where(detected_ids==0)) # Remove background id
    
    if len(detected_ids) == 0:
        return 0
    
    # Create consistent detected map
    consistent_detected_map = get_consistent_regions(detection_map, gt)    
    consistent_detected_map[consistent_detected_map>0] = 1
    # Create inconsistent detected map by removing the consistent region from the detection_map
    inconsistent_detected_map = detection_map.copy()
    inconsistent_detected_map[consistent_detected_map>0] = 0
    inconsistent_detected_map[inconsistent_detected_map>0] = 1

    # Calculates F1 score
    TP = np.sum(consistent_detected_map)
    FP = np.sum(inconsistent_detected_map)

    TPR = TP / (TP + FP )
    
    return TPR

def get_connected_components(cc_map,region_id):
    regions = np.zeros_like(cc_map)
    regions[cc_map == region_id] = 1
    regions_ids,ids = ndimage.label(regions)
    
    return regions_ids,ids

def get_consistent_regions(det_map,gt):
    """
    Return the consistent region map from detection map
    """
    
    gt_ids = np.unique(gt)
    gt_ids = np.delete(gt_ids,np.where(gt_ids==0)) # Remove background
    
    det_map_ids = np.unique(det_map)
    det_map_ids = np.delete(det_map_ids,np.where(det_map_ids==0)) # Remove background
    
    CTP = np.zeros_like(gt) # Consistent map
    
    # Loop troughtout the copied objects ids from the ground truth  
    for gid in gt_ids:
        gt_map, number_of_components = get_connected_components(gt,gid) # Get components with id = gid
        # from now on gt_map will locate each component relative from the forgery with id = gid
        # Each component is labeled with a different number
        
        # Now we need to check if the detection map has a region
        # that intersect more than one gt_map component
        
        for region_id in det_map_ids:
            region_map = np.zeros_like(det_map)
            region_map[det_map == region_id] = 1
            
            # get detected regions that intersection gt map
            det_inter_gt = np.logical_and(region_map, gt_map.copy())
            result = det_inter_gt * gt_map
            
            # if number of regions from results are greater or equal than 2, there was a match
            result_regions = np.unique(result)
            result_regions = np.delete(result_regions,0) # Remove background region
            if len(result_regions) >=2:
                result[result>0] = gid
                CTP += result.astype("uint8")
            
    return CTP

def is_object_consistent_with_gt(detected, groundtruth):
    """
    Check if the detected region is refering to the ground-truth
    Thus, we check if the connected component from the same label of a ground-truth region are 
    detected on the detected map.
    Since, some regions could have more than 2 connected componnets, or even just one,
    we considerate a dectection map valid as valid, if it can map at least two connected components ( or in the case 
    the gt has just one, the entire gt region for that object ID)
    """
    
    activation_map = detected.copy() # Activate map correspond to the regions activated by the solution
    gt = groundtruth.copy()
    
    gt[gt>0] = 1
    
    activation_map[activation_map > 1] = 1
   
    # for each region ( labeled with an ID) from the groundtruth we 
    # will check if there is a region on the activation map that correspond
    # to that region ( with an ID number that don't really need to match with the gt )
    gid = 1
    act_id = 1
    
    
    # Get connected compontents with the ID gid        
    cc_regions , cc_ids = get_connected_components(gt, gid)
    cc_regions = cc_regions.astype("uint8")


    # Initialize all connected componets as not Activated
    components_activated = [False for c in range(cc_ids)]

    # Loop on each conncected componnent from the groundtruth with id = gid
    for cc_id in range(1,cc_ids+1):
        
        gt_cc_map = np.zeros_like(cc_regions)
        gt_cc_map[cc_regions == cc_id] =1
        # Check if there is any intersection between the activation region act_id and the connecte component
        if  cv2.bitwise_and( activation_map, gt_cc_map ).sum() >0:
            components_activated[cc_id-1] = True

            
    # If there were at least 2 connected compontents from the ground-truth detected in the activatation map
    # or there is just one connected componet in the ground-truth and it was activated 
    # we return True
    if sum(components_activated) >= 2 or sum(components_activated) == len(components_activated):
        return True

    return False

def is_activation_map_consistent_with_gt(detected, groundtruth):
    """
    Check if the detected region is refering to the ground-truth
    Thus, we check if the connected component from the same label of a ground-truth region are 
    detected on the detected map.
    Since, some regions could have more than 2 connected componnets, or even just one,
    we considerate a dectection map valid as valid, if it can map at least two connected components ( or in the case 
    the gt has just one, the entire gt region for that object ID)
    """
    
    activation_map = detected.copy() # Activate map correspond to the regions activated by the solution
    gt = groundtruth.copy()
    
    gt_ids = np.unique(gt)
    gt_ids = gt_ids[gt_ids!=0] # Remove background id
   
    act_ids = np.unique(activation_map)
    act_ids = act_ids[act_ids!=0] # Remove background id
   
    # for each region ( labeled with an ID) from the groundtruth we 
    # will check if there is a region on the activation map that correspond
    # to that region ( with an ID number that don't really need to match with the gt )
    for gid in gt_ids:
        
        # Get connected compontents with the ID gid        
        cc_regions , cc_ids = get_connected_components(gt, gid)
        cc_regions = cc_regions.astype("uint8")
        
        
        # Loop on each activated region id
        for act_id in act_ids:
            
            # Initialize all connected componets as not Activated
            components_activated = [False for c in range(cc_ids)]
            
            # Loop on each conncected componnent from the region with id gid
            for cc_id in range(1,cc_ids+1):

                # Check if there is any intersection with the activation region act_id and the connected component
                act_map = np.zeros_like(activation_map)
                act_map[activation_map == act_id] = 1
                gt_cc_map = np.zeros_like(cc_regions)
                gt_cc_map[cc_regions == cc_id] =1
                
                
                """
                Debug
                plt.figure(figsize=(20,20))
                plt.subplot(1,3,1); plt.imshow(act_map)
                plt.subplot(1,3,2); plt.imshow(gt_cc_map)
                plt.subplot(1,3,3); plt.imshow(cv2.bitwise_and( act_map, gt_cc_map ))
                plt.show()
                """
                
                if  cv2.bitwise_and( act_map, gt_cc_map ).sum() >0:
                    components_activated[cc_id-1] = True
                    
            # If there were at least 2 connected compontents from the ground-truth detected in the activatation map
            # or there is just one connected componet in the ground-truth and it was activated 
            # we return True
            if sum(components_activated) >= 2 or sum(components_activated) == len(components_activated):
                return True
            
    return False
