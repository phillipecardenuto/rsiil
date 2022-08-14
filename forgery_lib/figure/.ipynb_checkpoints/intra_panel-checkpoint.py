"""
This module implements Intra Panel forgeries
Intra Panel forgeries are forgeries that occurs inside a unique panel of a Compound Figure
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
from .simple_figure import SimpleFigureForgery
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

class IntraPanelForgery(CompoundFigure):
    """
    This class inherit all methods from CompoundFigure to create a compound figure with
    Intra-Panel forgery ( just one panel of the figure is envolved in the forgery).
    You can use this class to insert a SinglePanelForgery into a panel of a Compound Figure.

    """




    def __init__(self, source_dataset, forgery_info, query, figure_template=None, template_dataset=None):
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
        drop_content = source_dataset[source_dataset['dataPath'] == query['dataPath']].index
        source_dataset = source_dataset.drop(drop_content)

        # Find valid Template from template dataset, if no specific template was input
        if figure_template is None:
            candidate_templates = self.find_valid_template(query, source_dataset, template_dataset)
            if candidate_templates != []:
                self.figure_template = random.choice(candidate_templates)
            else:
                raise IOError("No valid template was passed for the query Image")
        else:
            self.figure_template = figure_template

        # Validade the user input
        self._assert_input(source_dataset, forgery_info, self.figure_template,
                           query['dataPath'], query['dataGTPath'])

        self.forgery_image_path, self.forgery_obj_map_path = query['dataPath'], query['dataGTPath']

        # Load Template infos
        self.load_template_info(self.figure_template)

        # Load the source image dataset
        self._dataset = source_dataset

        # Forgery Class
        self.forgery_info = forgery_info
        self.forgery_class = forgery_info['forgery_class']

        if not self.forgery_class  in IntraPanelForgery.IMPLEMENTED_CLASSES:
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
            ar = self.pristine_image.shape[1] / self.pristine_image.shape[0]
            panel_id = self._select_panel_id(self.forgery_class, ar)
            self.panels[panel_id]['associated_img'] = self.pristine_image
            self.panels[panel_id]['class'] = '__PRISTINE__'
            self.panels[panel_id]['image_id'] = self.forgery_image_path


            # Assembly the artificial Figure using all panels
            super()._assembly_figure('pristine-panel')

        else:
            # Execute forgery on from one panel to another
            self._execute_intra_panel_forgery(forgery_info.copy(), self.forgery_class)

            # Assembly the artificial Figure using all panels
            super()._assembly_figure('intra-panel')

            self.gt_after_forgery = self.figure_groundtruth_forgery
            self.gt_before_forgery = self.figure_groundtruth_pristine

    def _assert_input(self, source_dataset, forgery_info, figure_template,
                      forgery_image_path, forgery_obj_map_path):
        """
        Check if the user's input is valid.
        """

        assert type(forgery_info) is dict, "forgery_info must be a dictonary"
        assert forgery_info['forgery_class'] in IntraPanelForgery.IMPLEMENTED_CLASSES
        assert 'function_name' in list(forgery_info.keys()), "'function_name' must be an arg forgery_info"
        if forgery_info['function_name'] != 'pristine':
            assert forgery_info['function_name'] in SimpleFigureForgery.IMPLEMENTED_FORGERIES ,  "Not implemented function"
            assert 'args' in list(forgery_info.keys()), "'args' must be an arg forgery_info "
            assert type(forgery_info['args']) is dict

        assert isinstance(source_dataset, pd.DataFrame), "source_dataset must be a dictionary"
        for atribute in ['class','dataGTPath', 'dataPath', 'AR']:
            assert atribute in source_dataset.columns, "source dataset must include 'class','dataGTPath', 'dataPath', 'AR' "
        for index, row in source_dataset.iterrows():
                dataPath = row['dataPath']
                assert os.path.isfile(dataPath) , f"{dataPath} is not a valid file"

        if not forgery_obj_map_path is None:
            assert  forgery_image_path, "forgery_obj_map_path, cannot be declared without a forgery_image_path"

        if not forgery_image_path is None:
            if forgery_info['function_name'] in SimpleFigureForgery.OBJ_MAP_NEEDED:
                assert forgery_obj_map_path, f"{forgery_info['function_name']} requires an object map as input"


    def find_valid_template(self, query, dataset, templates):

        # Since we are changing the dataset, we make a copy to avoid losing data
        aux_dataset = dataset.copy()
        return self._find_valid_template(query, aux_dataset, templates, figure_type='intra-panel')

    def get_forgery_metadata(self):
        """
        Return the metadata associated with the forgery.
        Since we are build a inter-panel figure,
        modality = get_forgery_modality(),
        figure_type = 'intra-panel'
        """
        return self._get_forgery_metadata(modality=self.get_forgery_modality(), figure_type='intra-panel')


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




