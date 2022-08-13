"""
This module implements InterPanel Forgeries
Inter Panel forgeries are forgeries that occurs within two or more panels
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
import pandas as pd
import json
import io
import string
import random
import os
from .compound_figure import CompoundFigure
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

class InterPanelForgery(CompoundFigure):
    """
    This class inherit all methods from CompoundFigure to create a compound figure with
    Inter-Panel forgery
    You can use this class to clone panels (with or without transformation), duplicate
    regions (with overlap forgery), or apply splicing_forgery in panels of a Compound Figure.
    
    """


 

    def __init__(self, source_dataset, forgery_info,
                 query=None, query_donor=None, query_host=None,
                 figure_template=None, template_dataset=None):
        """
        Parameters
        ----
        source_dataset: pandas object 
        Pandas dataframe containg all candidates images to be placed in the final Figure
        Columns of the pandas dataframe:
                ['class', 'DatasetRef', 'License', 'Link', 'ObjectMaskPath', 
                'RealDataPath', 'subset_tag', 'dataGTPath', 'dataPath', 'AR']
        
        forgery_info: <dict> 
        Dictonary with all info related with the type of forgery that the user wants to
        applied on a fake scietific figure.
            {
                function_name: <str>
                args:{
                     }
                multi_panel: <bool>  # Indicate if the foregery uses only one panel
                                   # or at least two panels
                forgery_class: <str>
            } 
        query: pandas dataframe 
            A row of the dataset containg the same columns of the source_dataset
            Use this input when the forgery function is not a splicing

         query_donor: pandas dataframe (for splicing forgeries)
            A row of the dataset containg the same columns of the source_dataset related to 
            the donor image from a forgery splicing.
            Whenever applying a splicing fucntion, make sure to explicit input this parameter
            InterPanel(source_dataset, forgery_info, template_dataset, query_donor=query_donor, query_host=query_host)

        query_host: pandas dataframe (for splicing forgeries)
            A row of the dataset containg the same columns of the source_dataset related to 
            the host image from a forgery splicing
            Whenever applying a splicing fucntion, make sure to explicit input this parameter
            InterPanel(source_dataset, forgery_info, template_dataset, query_donor=query_donor, query_host=query_host)
            
        figure_template: <string>
            Path to a Json file containing information relative to 
            the figure generation:
            Example:
            {

                panel_id: { bbox:{x0: <int> ,y0: <int>, x1: <int> ,y1: <int> },
                            class: <str>
                        },
                heighy: <int>,
                width: <int>
            }
        template_dataset: list(<string>)
            List of templates file path. This argument should be pass, if no specific template
            was already input, using the figure_template argument
            Each Template is designed to pointing the locations of its elements/panels in
            a figure
            {

                panel_id: { bbox:{x0: <int> ,y0: <int>, x1: <int> ,y1: <int> },
                            class: <str>
                        },
                heighy: <int>,
                width: <int>
            }
        """

        # Remove query from source_dataset
        if forgery_info['function_name'] == 'splicing_forgery':
            drop_content = source_dataset.loc[(source_dataset['dataPath'] == query_donor['dataPath']) |
                                              (source_dataset['dataPath'] == query_host['dataPath'])].index
        else:
            drop_content = source_dataset[source_dataset['dataPath'] == query['dataPath']].index
            
        source_dataset = source_dataset.drop(drop_content)
        
        # Forgery Infos
        self.forgery_info = forgery_info

        # Find valid Template from template dataset, if no specific template was input
        if figure_template is None:
            candidate_templates = self.find_valid_template(query, query_donor, query_host, source_dataset, template_dataset)
            if candidate_templates != []:
                self.figure_template = random.choice(candidate_templates)
            else:
                raise IOError("No valid template was passed for the query Image")
        else:
            self.figure_template = figure_template

        # Validade the user input
        self._assert_input(source_dataset, query, query_donor, query_host)
        
        # Load Template infos
        self.load_template_info(self.figure_template)

        # Load the source image dataset 
        self._dataset = source_dataset

        # Forgery Class
        self.forgery_class = forgery_info['forgery_class']
        
        if not self.forgery_class  in InterPanelForgery.IMPLEMENTED_CLASSES:
            raise NotImplementedError
            
        if forgery_info['function_name'] == 'pristine':
            # If forgery_class is pristine, we don't apply any forgery into the data
            # Therefore we just assert that there is a valid panel to the input image
            # from the forgery_image_path, then we assembly the final figure
            
            if self.forgery_image_path:
                # Check if there is a panel that could fit 
                # the image aspect ratio
                self._assert_panels_to_forgery(self.forgery_class,1)
            else:       
                # Select a image that could fit a random panel from forgery class
                self._select_forgery_image(self.forgery_class,1)
            # Insert Pristin image in a panel
            self.pristine_image_path = self.forgery_image_path
            self.pristine_image, _ = self.load_image(self.forgery_image_path)
            ar = self.pristine_image.shape[0] / self.pristine_image.shape[1]
            panel_id = self._select_panel_id(self.forgery_class, ar)
            self.panels[panel_id]['associated_img'] = self.pristine_image
            self.panels[panel_id]['class'] = '__PRISTINE__'
            self.panels[panel_id]['image_id'] = self.forgery_image_path
            
            
            # Assembly the artificial Figure using all panels
            self.__assembly_figure('pristine-panel')
        
        else:
            # Execute forgery on from one panel to another 
            self._execute_inter_panel_forgery(forgery_info.copy(), self.forgery_class)

            # Assembly the artificial Figure using all panels
            super()._assembly_figure('inter-panel')
            
            self.gt_after_forgery = self.figure_groundtruth_forgery
            self.gt_before_forgery = self.figure_groundtruth_pristine
        

    def _assert_input(self, source_dataset, query, query_donor, query_host):
        """
        Check if the user's input is valid. 
        """

        assert type(self.forgery_info) is dict, "forgery_info must be a dictonary"
        assert self.forgery_info['forgery_class'] in InterPanelForgery.IMPLEMENTED_CLASSES
        assert 'function_name' in list(self.forgery_info.keys()), "'function_name' must be an arg forgery_info"
        assert self.forgery_info['function_name'] in InterPanelForgery.IMPLEMENTED_FORGERIES ,  "Not implemented function"
        if self.forgery_info['function_name'] != 'pristine':
            assert 'args' in list(self.forgery_info.keys()), "'args' must be an arg forgery_info "
            assert type(self.forgery_info['args']) is dict
        
        assert isinstance(source_dataset, pd.DataFrame), "source_dataset must be a dictionary"
        for atribute in ['class','dataGTPath', 'dataPath', 'AR']:
            assert atribute in source_dataset.columns, "source dataset must include 'class','dataGTPath', 'dataPath', 'AR' "
        for index, row in source_dataset.iterrows(): 
                dataPath = row['dataPath']
                assert os.path.isfile(dataPath) , f"{dataPath} is not a valid file"    
        
        if self.forgery_info['function_name'] == 'splicing_forgery':
            self.forgery_info['donor_img_path'] = query_donor['dataPath'] 
            self.forgery_info['donor_obj_path'] = query_donor['dataGTPath'] 
            self.forgery_info['host_img_path'] = query_host['dataPath'] 
            self.forgery_info['host_obj_path'] = query_host['dataGTPath'] 
            assert os.path.isfile(query_host['dataPath']) , f"{query_host['dataPath']} is not a valid file"    
            assert os.path.isfile(query_donor['dataPath']) , f"{query_donor['dataPath']} is not a valid file"
            assert os.path.isfile(query_host['dataGTPath']) , f"{query_host['dataGTPath']} is not a valid file"    
            assert os.path.isfile(query_donor['dataGTPath']) , f"{query_donor['dataGTPath']} is not a valid file"        

        else:
            self.forgery_image_path, self.forgery_obj_map_path = query['dataPath'], query['dataGTPath']

            if not query['dataPath'] is None:
                if self.forgery_info['function_name'] in InterPanelForgery.OBJ_MAP_NEEDED:
                    assert query['dataGTPath'], f"{self.forgery_info['function_name']} requires an object map as input"

    
    def find_valid_template(self, query, query_donor, query_host, dataset, templates):
        
        # Since we are changing the dataset, we make a copy to avoid losing data
        aux_dataset = dataset.copy()
        if self.forgery_info['function_name'] == 'splicing_forgery':
            return self._find_splicing_template(query_donor, query_host, aux_dataset, templates)
        else:
            return self._find_valid_template(query, aux_dataset, templates, figure_type='inter-panel')

    def get_forgery_metadata(self):
        """
        Return the metadata associated with the forgery.
        Since we are build a inter-panel figure,
        modality = duplication,
        figure_type = 'inter-panel' 
        """
        return self._get_forgery_metadata(modality='duplication', figure_type='inter-panel')



    def _find_splicing_template(self, query_donor, query_host, dataset, templates):

        # Check the templates that has more than one panel with the same AR as the input image
        valid_templates = []
        for template_path in templates:
            with open(template_path,'r') as jf:
                template = json.load(jf)
            
            found = [0,0] # find template that matches for the donor (first position) and the host (last position)
            # check if we have data enough in the dataset to fill all panels
            image_allocated_to_figure = []
            for key, panel in template.items():
                if isinstance(panel,dict):
                    
                    panel['bbox']['y1'] = y1  = np.round(panel['bbox']['y1'] * template['orig_height']).astype(int)
                    panel['bbox']['y0'] = y0  = np.round(panel['bbox']['y0'] * template['orig_height']).astype(int)
                    panel['bbox']['x1'] = x1  = np.round(panel['bbox']['x1'] * template['orig_width']).astype(int)
                    panel['bbox']['x0'] = x0  = np.round(panel['bbox']['x0'] * template['orig_width']).astype(int)
                    
                    ar_panel = (y1-y0)/(x1-x0)
                    if (panel['class'] == query_donor['class']) and abs(ar_panel - query_donor['AR']) <= CompoundFigure.AR_TOLERENCE[query_donor['class']]:
                        found[0] +=1                                           
                    if (panel['class'] == query_host['class']) and abs(ar_panel - query_host['AR']) <= CompoundFigure.AR_TOLERENCE[query_host['class']]:
                        found[1] +=1                   

                if found[0] >=1 and found[1]>=1:
                    break
            if found[0] >=1 and found[1]>=1:
                valid_templates.append(template_path)
                
    
    # Check if we have enought data to fill the template
        invalid_templates = []
        for template_path in valid_templates:
            # pandas dataset
            candidate_dataset = dataset.copy()
            with open(template_path,'r') as jf:
                template = json.load(jf)

            # Create a pandas dataframe to optimize the search of images in the dataset
            candidate_panels = pd.DataFrame(columns=['class','AR','ID','filled']) 

            for key, panel in template.items():
                if isinstance(panel,dict):
                    panel['bbox']['y1'] = y1  = np.round(panel['bbox']['y1'] * template['orig_height']).astype(int)
                    panel['bbox']['y0'] = y0  = np.round(panel['bbox']['y0'] * template['orig_height']).astype(int)
                    panel['bbox']['x1'] = x1  = np.round(panel['bbox']['x1'] * template['orig_width']).astype(int)
                    panel['bbox']['x0'] = x0  = np.round(panel['bbox']['x0'] * template['orig_width']).astype(int)
                    ar_panel = (y1-y0)/(x1-x0)

                    filled = False
                    if panel['class'] in ['Others', 'Graphs']:
                        filled = True
                    candidate_panels = candidate_panels.append({"class":panel['class'],
                                            "AR":ar_panel,
                                            "ID":key,
                                            "filled":filled}, ignore_index=True)

            # Reserve one panel to the donor and other to the host
            index = candidate_panels.loc[(candidate_panels['class']==query_donor['class']) & 
                                                ((query_donor['AR']-CompoundFigure.AR_TOLERENCE[query_donor['class']]) <= candidate_panels.AR) &
                                                (candidate_panels.AR <= (query_donor['AR']+CompoundFigure.AR_TOLERENCE[query_donor['class']]))].iloc[0:1].index
            candidate_panels.loc[index,'filled'] = 'query'
            candidate_panels.loc[index,'class'] = 'query'
            index = candidate_panels.loc[(candidate_panels['class']==query_host['class']) & 
                                                ((query_host['AR']-CompoundFigure.AR_TOLERENCE[query_host['class']]) <= candidate_panels.AR) &
                                                (candidate_panels.AR <= (query_host['AR']+CompoundFigure.AR_TOLERENCE[query_host['class']]))].iloc[0:1].index
            candidate_panels.loc[index,'filled'] = 'query'
            candidate_panels.loc[index,'class'] = 'query'

            # Try to fill the rest of the panels
            for index, row_panel in candidate_panels.loc[candidate_panels['filled'] == False].iterrows():
                cand_img = candidate_dataset.loc[(candidate_dataset['class']==row_panel['class']) & ((candidate_dataset.AR-CompoundFigure.AR_TOLERENCE[row_panel['class']]) <=row_panel['AR']) &
                                                (row_panel['AR'] <= candidate_dataset.AR+CompoundFigure.AR_TOLERENCE[row_panel['class']])].iloc[0:1]
                if len(cand_img) == 0:
                    invalid_templates.append(template_path)
                    break
                else:
                    candidate_dataset = candidate_dataset.drop(cand_img.index)

        for template_path in invalid_templates:
            valid_templates.remove(template_path)
                    
        return valid_templates


        