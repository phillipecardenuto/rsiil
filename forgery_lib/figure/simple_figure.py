"""
Simple Forgery Figures:
Creates a Simple Forgery figures from a implemented forgery 

-------------------------------------
Author: Phillipe Cardenuto
email: phillipe.cardenuto@gmail.com
--------------------------------------
Version:                           1.1
Version Date:               Nov , 2020

This is file is part of Forgery Library
"""

import numpy as np
import cv2
import itertools
from PIL import Image, ImageDraw, ImageFont, ImageOps
import matplotlib.pyplot as plt
import json
import io
import string
import random
import os
from ..retouching import (retouching_contrast,
                         retouching_brightness,
                         retouching_brightness_contrast,
                         retouching_blurring)
from ..cleaning import(  cleaning_using_bf,
                         cleaning_with_inpainting)
from ..duplication import (copy_move_forgery,
                        random_copy_move,
                        splicing_forgery,
                        overlap_forgery,
                        simple_copy)

class SimpleFigureForgery():
    """
    This class implements Simple Figure Forgery.
    To organize the forgeries implemented in the library, we insert all forgery in 
    the IMPLEMENTED_FORGERIES list. If you apply any third-party forgery function, 
    remind to include the name of this function in this list.
    
    """
    IMPLEMENTED_FORGERIES = [ 'retouching_contrast', 'retouching_brightness',
                         'retouching_brightness_contrast', 'retouching_blurring',
                         'cleaning_using_bf', 'cleaning_with_inpainting',
                         'copy_move_forgery','random_copy_move', 'splicing_forgery',
                         'overlap_forgery', 'simple_copy'
                         ]
    OBJECT_GROUNDTRUTH = ['cleaning_using_bf','copy_move_forgery',
                           'random_copy_move']

    MULTIPLE_GROUNDTRUTH = ['splicing_forgery','overlap_forgery', 'simple_copy']

    DUPLICATION_MODALITY = ['copy_move_forgery','random_copy_move',
                            'splicing_forgery','overlap_forgery', 'simple_copy'] 

    RETOUCHING_MODALITY = [ 'retouching_contrast', 'retouching_brightness',
                         'retouching_brightness_contrast', 'retouching_blurring']

    CLEANING_MODALITY = [ 'cleaning_using_bf', 'cleaning_with_inpainting']

    OBJ_MAP_NEEDED = ['retouching_contrast', 'retouching_brightness',
                         'retouching_brightness_contrast', 'retouching_blurring',
                         'cleaning_using_bf', 'cleaning_with_inpainting',
                         'copy_move_forgery','random_copy_move']
    def __init__(self, img_path, forgery_info, objs_maks_path=''):
        """
        Parameters
        ----
        img_path: <str>
            Path to the source image that will suffer forgery
            
        forgery_info: <dict> 
            This dictonary has intructions of how to apply the forgery
            {
                function_name: str
                args:{ ... }
            } 

            In case of applying splicing forgery, use the forgery_info to 
            pass all arguments related to the host/donor images and their ground-truth,
            and pass img_path=None:

            {
                donor_img_path : str,
                donor_obj_path:  str,
                host_img_path: str,
                host_obj_path: str
            }


        objs_maks_path: (opitional) <str>
            Path to the groundtruth object mask
        """
        ###########
        # Asserts #
        ###########
        if forgery_info['function_name'] == 'splicing_forgery':
            assert 'donor_img_path' in list(forgery_info.keys()), "splicing forgery function must have a donor_img_path key"
            assert 'donor_obj_path' in list(forgery_info.keys()), "splicing forgery function must have a donor_obj_path key"
            assert 'host_img_path' in list(forgery_info.keys()), "splicing forgery function must have a host_img_path key"
            assert 'host_obj_path' in list(forgery_info.keys()), "splicing forgery function must have a host_obj_path key"
            assert img_path is None, "splicing forgery function must have img_path = None"

            donor_img_path, donor_obj_path =  forgery_info['donor_img_path'], forgery_info['donor_obj_path']
            host_img_path, host_obj_path  = forgery_info['host_img_path'], forgery_info['host_obj_path'] 

            self.donor_img, self.donor_obj_map = self.__load_image(donor_img_path, donor_obj_path)
            self.host_img, self.host_obj_map = self.__load_image( host_img_path, host_obj_path)

        else:
            assert os.path.isfile(img_path),  f"{img_path} is not a file"
            assert type(forgery_info) is dict, "forgery_info must be a dictonary"
            assert 'function_name' in list(forgery_info.keys())
            assert 'args' in list(forgery_info.keys())
            assert type(forgery_info['args']) is dict
            assert forgery_info['function_name'] in SimpleFigureForgery.IMPLEMENTED_FORGERIES
            if forgery_info['function_name'] in SimpleFigureForgery.OBJ_MAP_NEEDED:
                assert objs_maks_path != ''

            self.img_path = img_path
            self.objs_mask_path = objs_maks_path
            
            self.src_image, self.objs_mask = self.__load_image(img_path,objs_maks_path )


        if not 'forgery_class' in forgery_info.keys():
            forgery_info['forgery_class'] = ""
        
        self.forgery_info = forgery_info 

        #################
        # Apply Forgery #
        #################
        self.__apply_forgery()
        
    def __apply_forgery(self):
        """
        Using the forgery type apply the forgery on the pristine image
        """
        function_name = self.forgery_info['function_name']
        args = self.forgery_info['args']

        if function_name == 'splicing_forgery':
            args['donor'] = self.donor_img.copy()
            args['donor_map'] = self.donor_obj_map.copy()
            args['host'] = self.host_img.copy()
            args['host_map'] = self.host_obj_map.copy()

        else:
            args['img'] = self.src_image.copy()
            if self.forgery_info['function_name'] in SimpleFigureForgery.OBJ_MAP_NEEDED:
                args['objs_map'] = self.objs_mask.copy()
        
        # Abstract forgery call
        forgery_call = f"{function_name}(**args)"
        output = eval(forgery_call)
        self._organize_groudtruth(output)
        
       

    def __load_image(self, img_path, objs_mask_path=''):
        """
        Load an image from a path
        TODO: Deal with extensions problems that PIL can't
        """
        img = np.array(Image.open(img_path))
        if objs_mask_path:
            objs_map = np.array(Image.open(objs_mask_path).convert('L'))
        else:
            objs_map = None
        return img, objs_map


    def _organize_groudtruth(self, output):
        function_name = self.forgery_info['function_name']

        if function_name in SimpleFigureForgery.OBJECT_GROUNDTRUTH:
            
            f_img, f_gt, src_gt = output
            if f_img is None:
                raise RuntimeError(f"The mehtod{function_name} not generated a forgery image.")
            self.forgery_img = f_img
            self.gt_after_forgery = f_gt
            self.gt_before_forgery = src_gt

        elif function_name == 'simple_copy':
            src_img, src_gt, f_img, f_gt = output

            if f_img is None:
                raise RuntimeError(f"The mehtod{function_name} not generated a forgery image.")
            self.forgery_img = f_img
            self.gt_after_forgery = f_gt
            self.gt_before_forgery = src_gt

        elif function_name == 'overlap_forgery':
            region1, region1_gt, region2, region2_gt, = output

            if region1 is None or region2 is None:
                raise RuntimeError(f"The mehtod{function_name} not generated a forgery image.")
            self.region1 = region1
            self.region1_gt = region1_gt
            self.region2 = region2
            self.region2_gt = region2_gt

        elif function_name == 'splicing_forgery':
            f_img, donor_gt_map, host_gt_map = output

            if f_img is None:
                raise RuntimeError(f"The mehtod{function_name} not generated a forgery image.")
            self.forgery_img = f_img
            self.donor_objects_gt = donor_gt_map
            self.host_objects_gt = host_gt_map

        else:
            # Neither Objects or regions change its location before or after the forgery
            f_img, f_gt = output
            if f_img is None:
                raise RuntimeError(f"The mehtod{function_name} not generated a forgery image.")
            self.forgery_img = f_img
            self.gt_after_forgery = f_gt
            self.gt_before_forgery = f_gt

    def get_forgery_metadata(self):
        """
        Return the metada involved in the forgery 
        """
        metadata = {}
        forgery_info_args = self.forgery_info['args'].copy()
        

        if self.forgery_info['function_name'] == 'splicing_forgery':
            del forgery_info_args['donor']
            del forgery_info_args['donor_map'] 
            del forgery_info_args['host']
            del forgery_info_args['host_map']
        else:
            del forgery_info_args['img']
            if self.forgery_info['function_name'] in SimpleFigureForgery.OBJ_MAP_NEEDED:
                del forgery_info_args['objs_map']

        metadata['forgery_info'] = {
                'class': self.forgery_info['forgery_class'],
                'function_name': self.forgery_info['function_name'],
                'args': forgery_info_args,
                "figure_type": "simple",
                'modality': self.get_forgery_modality()
        }
        
        if self.forgery_info['function_name'] == 'splicing_forgery':
            metadata["figure_annotations"] = { "donor_id": self.forgery_info['donor_img_path'],
                                            "host_id":  self.forgery_info['host_img_path'],
                                                "height": self.host_img.shape[0],
                                                "width": self.host_img.shape[1],
                                                "image_id": self.forgery_info['donor_img_path']
            }       

        else:
            metadata["figure_annotations"] = {
                                "height": self.src_image.shape[0],
                                "width": self.src_image.shape[1],
                                "image_id": self.img_path
                    }

        return metadata
    

    def get_forgery_modality(self):
        """
        Return the type of the forgery based on its name
        """
        if self.forgery_info['function_name'] in SimpleFigureForgery.DUPLICATION_MODALITY:
            return 'duplication'
        elif self.forgery_info['function_name'] in SimpleFigureForgery.RETOUCHING_MODALITY:
            return 'retouching'
        elif self.forgery_info['function_name'] in SimpleFigureForgery.CLEANING_MODALITY:
            return 'cleaning'
        else:
            return 'Not Implemented'
