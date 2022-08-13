"""
This module implements Compound figures Forgeries
Compound Figures:
    Typically scientific figures are compound by panels. If the user choose this level,
    a template figure with the design of the figure is required. Also, the images that
    would be included in the panel are required.

    Using a JSON file as template, this class creates
    a scientific figure applying common forgery saw on retracted papers.

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
import copy
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

class CompoundFigure():

    """
    This is an abstract class to create Compound figures
    """


    VERBOSITY_LEVEL3_CLASSES = ['Microscopy']

    # Aspect Ratio Tolerance
    # Since for different classes the tolerance for fit in a panel can be different
    # We define the tolerence per class

    AR_TOLERENCE = {
            'WesternBlot': 2,
            'Microscopy': 4
    }

    OBJ_MAP_NEEDED = ['retouching_contrast', 'retouching_brightness',
                         'retouching_brightness_contrast', 'retouching_blurring',
                         'cleaning_using_bf', 'cleaning_with_inpainting',
                         'copy_move_forgery','random_copy_move']



    IMPLEMENTED_FORGERIES = [ 'retouching_contrast', 'retouching_brightness',
                         'retouching_brightness_contrast', 'retouching_blurring',
                         'cleaning_using_bf', 'cleaning_with_inpainting',
                         'copy_move_forgery','random_copy_move', 'splicing_forgery',
                         'overlap_forgery', 'simple_copy', 'pristine']

    # List of function that need at least two panels to be successfuly applied
    MULTI_PANEL_FORGERY = ['simple_copy','splicing_forgery','overlap_forgery']

    IMPLEMENTED_CLASSES = ['Microscopy', 'WesternBlot']

    def __init__(self):
        pass
        """
        This class is used to construct CompoundFigure
        as an abstract class

        self.template_image = None
        self.width = 0
        self.height = 0
        self.panels = None
        self.forgery_class = ""
        raise "Do not instance any object from this class"
        """

    def insert_panel_in_figure(self,figure,panel_figure, bbox):

        # Read image
        if type(panel_figure) is str:
            panel_img = np.array(Image.open(panel_figure))
        else:
            panel_img = panel_figure
        # If panel_img is a GrayScale image
        if len(panel_img.shape) < 3:
            panel_img = cv2.merge((panel_img,panel_img,panel_img))
        # If panel_img has alpha
        elif panel_img.shape[-1] > 3:
            panel_img = panel_img[:,:,:3]

        bbox_img = (slice(bbox['y0'],bbox['y1']),
                slice(bbox['x0'],bbox['x1']))
        height, width = figure[bbox_img].shape[:2]
        panel_img = cv2.resize(panel_img,(width,height))
        figure[bbox_img] = panel_img
        return figure

    def make_plt_image(self,fig,width, height):
        """
        Creates a image from a matplotlib graph
        """
        io_buf = io.BytesIO()
        plt.savefig(io_buf, format='png', bbox_inches='tight')
        io_buf.seek(0)
        img = Image.open(io_buf)
        # img_arr = np.reshape(np.frombuffer(io_buf.getvalue(), dtype=np.uint8),
                            #  newshape=(int(fig.bbox.bounds[3]), int(fig.bbox.bounds[2]), -1))
        img_arr = np.array(img)
        io_buf.close()
        return img_arr[:,:,:3]


#-----------------------------------------------------------------------------------------------------#
#                                                                                                     #
#                                 METHODS RELATED WITH FAKE GRAPH                                     #
#                                                                                                     #
#-----------------------------------------------------------------------------------------------------#


    def generate_fake_bar(self, figshape,dpi=500):
        """
        With some sort of randomness creates a bar graph
        REF: https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py
        """

        # width of the bars
        barWidth = 0.1
        # Get random number of bars and datas
        number_of_group_data = np.random.randint(2,5)
        number_of_data= np.random.randint(2,5)

        # Define the possible colors of the bars
        COLOR_LIST = ['black','dimgray','gray','darkgray','silver',
                                                   'lightgray','gainsboro','whitesmoke',
                                                   'slategrey','lightslategray','white']
        # Shuffle list randomly
        np.random.shuffle(COLOR_LIST)
        colors = itertools.cycle(COLOR_LIST)

        # Flag to decide if the bar will be colored or not
        COLOR = np.random.randint(2)

        # Define possible patters to b&W bars
        PATTERNS = [ '',"|" , "\\" , "/" , "+" , "-", ".", "*","x", "o", "O" ]
        np.random.shuffle(PATTERNS)
        PATTERN_LIST = itertools.cycle(PATTERNS)

        # Generate random bars data
        bars = []
        for j in range(number_of_data):
            d = [ [random.randrange(1,10,1) , random.randrange(1,100,1)/100 ]# Heigh and errors
                 for i in range(number_of_group_data)]
            bars.append(d)
        bars = np.array(bars)

        # Get position of each data
        x = np.arange((number_of_group_data))

        # Create high resolution graph
        h = (figshape[0])
        w = (figshape[1]*5)/h
        fig,ax = plt.subplots(figsize=(w, 5),dpi=dpi,
                              frameon=False)

        groups_rects = []

        # Generate Bar plot
        # For each set data bar, plot it at the
        # right position
        for index,bar in enumerate(bars):
            # Get height and error of the bar
            height = bar[:,0]
            yerr = bar[:,1]

            # Get X position of each bar
            pos_bar = x - (number_of_data*barWidth/2) + index*barWidth
            # Plot bar with color
            if COLOR:
                groups_rects.append(
                        ax.bar(pos_bar,height,
                        yerr=yerr,
                        width=barWidth,
                        color=next(colors),
                        edgecolor='black',
                       capsize=8, label=f'FkD{index+1}')
                        )
            # Plot bar with hatch
            else:
                hatch = next(PATTERN_LIST)
                hatch = 3*hatch
                groups_rects.append(
                        ax.bar(pos_bar,height,
                        yerr=yerr,
                        width=barWidth,
                        color='white',
                        edgecolor='black',
                        hatch=hatch,
                        capsize=8, label=f'FkD{index+1}')
                        )

        # insert x-ticks addressing different name
        # for each set of bars
        plt.xticks(x,
                   [f'Fake{chr(65+index)}' for index in x])

        plt.ylabel('FakeBar')
        ax.legend(loc='upper right')
        # plt.title('FakeData')

        # Make image from plot
        img = self.make_plt_image(fig,w,h)
        plt.close()

        return img

    def generate_fake_plot(self,figshape,dpi=500):
        """
        With some sort of randomness creates a line graph
        """
        # High resolution Graph
        h = (figshape[0])
        w = (figshape[1]*5)/h
        fig,ax = plt.subplots(figsize=(w, 5),dpi=dpi,
                              frameon=False)

        # Randomly select the number of lines from 1 to 5
        number_of_lines = np.random.randint(1,5)
        # Randomly select the number of points from 3 to 10
        number_of_points = np.random.randint(3,10)

        markers = itertools.cycle(np.random.choice(["v","^","<",">","s","p","X","D","o"],number_of_lines))
        linestyles =  itertools.cycle(np.random.choice(['-','--', '-.', ':',],number_of_lines))
        # Get Lines
        lines = np.random.rand(number_of_lines,number_of_points,2)
        color = np.random.choice([None,'black'])
        for id_, line in enumerate(lines,start=1):
            line[:,0] = sorted(line[:,0])
            x = line[:,0]

            y = sorted(line[:, 1],reverse=np.random.choice([0,1]))
            x = (x - min(x)) / (max(x) - min(x))

            ax.scatter(x,y,s=100,marker=next(markers),label=f"FakeLine {id_}",color=color)

            ax.plot(x,y,linestyle=next(linestyles),color=color)

        ax.legend(loc='upper right')

        # Make image from plot
        img = self.make_plt_image(fig,w,h)
        plt.close()
        return img

#-----------------------------------------------------------------------------------------------------#
#                                                                                                     #
#                           METHODS RELATED WITH INDICATIVE CAPTION                                   #
#                                                                                                     #
#-----------------------------------------------------------------------------------------------------#

    def apply_verbose_text(self, figure, verbosity_level,
                           font_name='forgery_lib/figure/fonts/times-new-roman-14.ttf'):
        """
        Apply random text in the figure with some level of verbosity.
        As the level of verbose is higher more occlusion --'informative text'--
        will be added to the figure.

        Level 1: Includes only indicative letters around each panel
        Level 2: Includes all features from Level 1 and random words
                 around each panel
        Level 3: Includes all features from Level 2 and letters inside
                 the panels. (Only Panels which class are inside the list CompoundFigures.VERBOSITY_LEVEL3_CLASSES
                 will receive LEVEL 3 of verbosity). This feature will avoid images such as western blots to receive
                 a letter inside their panels.
        Parameters:
        ----------
        figure: numpy image
        verbosity_level: int in range of [0,3]
        """
        # Reset all letters in template image
        self.template_image[self.template_image == 255] =0
        if verbosity_level > 0:
            figure = self._verbosity_level1(figure,font_size=30,font_name=font_name)
        if verbosity_level > 1:
            figure = self._verbosity_level2(figure,text_length=5,font_name=font_name,
                                                  font_size=30)
        if verbosity_level > 2:
            figure = self._verbosity_level3(figure,font_name=font_name)
        return figure



    def _verbosity_level3(self,
                              figure,
                              font_name="forgery_lib/figure/fonts/times-new-roman-14.ttf",
                              font_size=25):
        """
        Find each level3 panel of the figure, based on the template, and insert
        a indicative letter inside its panel, e.g. (A).
        """

        # Check whether the foregery panel allow to receive LEVEL3, otherwise will not include it to
        # the list of panels that would receive a letter inside the panel
        if self.forgery_class in CompoundFigure.VERBOSITY_LEVEL3_CLASSES:
            list_level3_classes =  CompoundFigure.VERBOSITY_LEVEL3_CLASSES + ["__FORGERY__","__PRISTINE__"]
        else:
            list_level3_classes = CompoundFigure.VERBOSITY_LEVEL3_CLASSES

        level3_panels = []
        for _, panel in self.panels.items():
            if panel['class'] in list_level3_classes :

                bbox = panel['bbox']
                level3_panels.append(bbox)

        # if none level3_panels were found
        if level3_panels == []:
            return figure

        # Sort lexicographic the coordinate of free_spaces
        level3_panels.sort(key=lambda bx: (bx['x0'],bx['y0']))
        alphabet = itertools.cycle(string.ascii_lowercase)

        for panel in level3_panels:
            bbox_text = (slice(panel['y0'],panel['y0']+font_size,None ),
                         slice(panel['x0'],panel['x0']+font_size,None))
            figure = self._insert_text(figure,
                            next(alphabet),
                            panel['x0'],
                            panel['y0'],
                            font_name=font_name,
                            font_size=font_size,
                            color_font=self._get_color_font(figure,bbox_text))
        return figure

    def _check_panel_legend(self, font_size, x0,y0):
        """
        This method is use during indicative text level 1, to insert a letter
        to indicate a panel.
        Check if the neighborhood of a panel (x0,y0) border has any overlap with
        other panel.
        Here we consider just the right side and the top side as neighborhood
        """
        # Check if there the free space apart of the x,y fall outside the image
        if x0-font_size -1< 0 or y0+font_size > self.template_image.shape[0]:
            return False

        # Check if the right bbox has overlap with other panel
        right_side_bbox = (slice(y0,y0+font_size,None),
                           slice(x0-font_size,x0,None))

        caption_area = self.template_image[right_side_bbox]
        if np.count_nonzero(caption_area) != 0 :
            return False
        else:
            self.template_image[right_side_bbox] = 255

        return True

    def _verbosity_level1(self,
                                 figure,
                                 font_name="forgery_lib/figure/fonts/times-new-roman-14.ttf",
                                 font_size=45):
        """
        Find each panel that its right upper space is free to insert
        a caption letter. After this, insert a lexicographic sorted caption on
        each free space.
        """

        free_spaces = []
        for _, panel in self.panels.items():

            bbox = panel['bbox']
            x0,y0 = bbox['x0'], bbox['y0']
            if self._check_panel_legend(font_size,x0,y0):
                free_spaces.append((x0-font_size,y0))

        # If there is none free_space for caption, return
        if len(free_spaces) == 0 :
            return figure

        # Sort lexicographic the coordinate of free_spaces
        free_spaces.sort(key=lambda n: n[::-1])
        # Using include the letter text fon the free_space
        # free_spaces [ (bbox, 'A'), ... ]
        free_spaces = [(f, chr(65+index)) for index, f in enumerate(free_spaces)]

        for f in free_spaces:
            figure = self._insert_text(figure=figure,
                                       text=f[1],
                                       x=f[0][0],
                                       y=f[0][1],
                                       font_name=font_name,
                                       font_size=font_size)
        return figure

    def _verbosity_level2(self,
                                   figure,
                                   text_length,
                                   font_name="forgery_lib/figure/fonts/times-new-roman-14.ttf",
                                   font_size=10):
        """
        Insert a random word description with len text_size in a free region (that do not overlaps other panels) around
        the panels.
        """

        # Find all free spaces
        free_spaces = []
        for _, panel in self.panels.items():

            bbox = panel['bbox']
            region = self._get_free_border_region(font_size,text_length,bbox)
            if region:
                free_spaces.append((bbox,region))

        # If there is none free_space for caption, return
        if len(free_spaces) == 0 :
            return figure

        for f in free_spaces:
            figure = self._insert_description_text(figure=figure,
                                               text_length=text_length,
                                               bbox=f[0],
                                               direction=f[1],
                                               font_name=font_name,
                                               font_size=font_size)
        return figure

    def _insert_description_text(self,figure,
                             text_length,
                             bbox,
                             direction,
                             font_name="forgery_lib/figure/fonts/times-new-roman-14.ttf",
                             font_size=20):
        """
        Insert random verbose text inside a bbox, rotating the bbox according
        the direction param
        # REF: https://stackoverflow.com/questions/245447/how-do-i-draw-text-at-an-angle-using-pythons-pil
        """

        fig_pil = Image.fromarray(figure)

        # Height and width of the text
        h = font_size
        w = font_size*(text_length//2+1)


        # create image with size w,h
        txt_pil = Image.new("L",(w,h),255)

        # Generate text
        text = self._get_random_text(text_length)

        # Create draw
        draw = ImageDraw.Draw(txt_pil)
        # Get Font object
        font = ImageFont.truetype(font_name, size=font_size)
        draw.text((0,0), text, font=font, fill='rgb(0,0,0)')

        if direction == 'left':
            paste_box = (bbox['x0']-font_size*(1+text_length//2)-1,
                         bbox['y0'])
            fig_pil.paste( txt_pil,paste_box )
            # Mark Region on template image
            paste_area = (slice(paste_box[1],paste_box[1]+txt_pil.height,None),
                          slice(paste_box[0],paste_box[0]+txt_pil.width,None))
            self.template_image[paste_area] = 255

        elif direction == 'right':
            rot=txt_pil.rotate(270,expand=1)
            paste_box = (bbox['x1']+1,
                         bbox['y0'])
            fig_pil.paste( rot, paste_box )
            # Mark Region on template image
            paste_area = (slice(paste_box[1],paste_box[1]+rot.height,None),
                          slice(paste_box[0],paste_box[0]+rot.width,None))
            self.template_image[paste_area]

        elif direction == 'top':
            paste_box= (bbox['x0'],
                        bbox['y0']-font_size-1)
            fig_pil.paste(txt_pil, paste_box)
            # Mark Region on template image
            paste_area = (slice(paste_box[1],paste_box[1]+txt_pil.height,None),
                          slice(paste_box[0],paste_box[0]+txt_pil.width,None))
            self.template_image[paste_area] = 255

        elif direction == 'bottom':
            paste_box = (bbox['x0'],
                         bbox['y1']+1)
            fig_pil.paste(txt_pil, paste_box)
            # Mark Region on template image
            paste_area = (slice(paste_box[1],paste_box[1]+txt_pil.height,None),
                          slice(paste_box[0],paste_box[0]+txt_pil.width,None))
            self.template_image[paste_area] = 255

        return np.array(fig_pil)

    def _insert_text(self,
                    figure,
                    text,
                    x,y,
                    font_name="forgery_lib/figure/fonts/times-new-roman-14.ttf",
                    font_size=45,
                    color_font='rgb(0,0,0)'):
        """
        Insert text in the figure in the location of bbox.
        This function is used only for caption, legend propouse

        Return
        Figure with the text inserted
        """


        fig = Image.fromarray(figure)
        draw = ImageDraw.Draw(fig)
        font = ImageFont.truetype(font_name, size=font_size)
        draw.text((x,y), text,
                        fill=color_font, font=font)

        return np.array(fig)

    def _get_free_border_region(self,font_size, text_length, bbox):
        """
        Get the border of the panel(bbox) that a caption with font_size and
        text_length could be display.
        The preference of diplay is top, left, right, bottom
        """
        x0,y0 = bbox['x0'],bbox['y0']
        x1,y1 = bbox['x1'],bbox['y1']

        # Check top border
        if y0 - font_size -1 >0 :
            check_area = ( slice(y0-font_size-1,y0,None),
                           slice(x0,x1,None))
            # Check if the area overlap with any panels region
            caption_area = self.template_image[check_area]
            if np.count_nonzero(caption_area) == 0 :
                return "top"

        # Check left border
        if x0 - font_size*text_length - 1 > 0:
            check_area = ( slice(y0,y1,None),
                           slice(x0- font_size*(1+text_length//2) -1, x0, None))
            # Check if the area overlap with any panels region
            caption_area = self.template_image[check_area]
            if np.count_nonzero(caption_area) == 0 :
                return "left"

        # Check right border
        if x1 + font_size + 1 < self.template_image.shape[1] :
            check_area = ( slice(y0,y1,None),
                           slice(x1, x1+ font_size+1, None))
            # Check if the area overlap with any panels region
            caption_area = self.template_image[check_area]
            if np.count_nonzero(caption_area) == 0 :
                return "right"

        # Check bottom border
        if y1 + font_size + 1 < self.template_image.shape[0]:
            check_area = ( slice(y1,y1+font_size+1,None),
                           slice(x0, x1, None))
            # Check if the area overlap with any panels region
            caption_area = self.template_image[check_area]
            if np.count_nonzero(caption_area) == 0 :
                return "bottom"


        return None

    def _get_random_text(self,size=4):
        """
        Using string library and numpy random
        generate a random string of length 'size'
        """
        return ''.join( np.random.choice(
                           list(string.ascii_letters+string.digits),
                           size=size
                        )
                      )

    def _get_color_font(self,figure,bbx):
        """
        Based on the histogram of the bbx, return a contrast color
        """

        figure_pil= Image.fromarray(figure)
        if not 'L' in figure_pil.getbands():
            figure_pil.convert("L")

        # Check histogram of bouding box
        img_array  = np.array(figure_pil)
        img_array = img_array[bbx]
        hist, bins = np.histogram(img_array)
        # If value of hist concentrate in a dark region
        if bins[np.argmax(hist)] < 100:
            return "rgb(255,255,255)"
        # If value of hist concentrate in bright region
        return "rgb(0,0,0)"


#-----------------------------------------------------------------------------------------------------#
#                                                                                                     #
#                           METHODS RELATED WITH FIGURE TEMPLATE                                      #
#                                                                                                     #
#-----------------------------------------------------------------------------------------------------#

    def _create_template_image(self):
        """
        Creates a template image aimming any sort of visulization or
        auxiliary usage during the figure creation
        """

        template_image = np.zeros((self.height,self.width))
        figures_classes = []

        for _, panel in self.panels.items():
            if panel['class'] in figures_classes:
                class_index = figures_classes.index(panel['class'])+1
            else:
                figures_classes.append(panel['class'])
                class_index = figures_classes.index(panel['class'])+1

            bbox = (slice(panel['bbox']['y0'],panel['bbox']['y1'],None),
                    slice(panel['bbox']['x0'],panel['bbox']['x1'],None))
            template_image[bbox] = class_index

        return template_image


    def validate_template(self):
        pass


    def _validate_number_of_panels(self, forgery_class, candidate_panels_ids, minum_number_of_panels=1):
        """
        Check whether the template layout has a minimum number of panels of class 'forgery_class'
        """

        # check if candidates panels are enough
        if len(candidate_panels_ids) < minum_number_of_panels:
            assert False, f"The Layout of the template does not have enough panels from class {forgery_class}\
                            to complete the forgery"

        # Check if there are n minum number of panels with the same aspect rate
        while len(candidate_panels_ids) >= minum_number_of_panels:

            cand_id = random.choice(candidate_panels_ids)
            candidate_panels_ids.remove(cand_id)
            cand_panel = self.panels[cand_id]

            accepted_panels = 1
            for _, panel in self.panels.items():
                if panel['class'] == forgery_class and  \
                   np.absolute( cand_panel['ar'] - panel['ar'] ) <= 2 and \
                    panel['associated_img'] is None:
                    accepted_panels+=1
                    if accepted_panels >= minum_number_of_panels:
                        return cand_id

        assert False, f"The Layout of the template does not have enough panels from class {forgery_class}\
                        to complete the forgery"


    def load_template_info(self,figure_template):
        """
        Load a Template Json File from the disk and create an empty image from it
        """

        if type(figure_template) is str:
            assert os.path.isfile(figure_template), "If figure_template is a string it must be a path to a json file"
            with open(figure_template, 'r') as rj:
                self.template = json.load(rj)
        else:
            # Make a deep copy of the figure template input
            import copy
            self.template = copy.deepcopy(figure_template)
        # Get the height of the final figure
        # Since we are preserving the aspect rate of the template
        # the width will be rescaled according the height users input
        if not 'height' in self.template.keys():
            self.height = int(self.template['orig_height'])
        else:
            self.height = int(self.template['height'])

        self.width = np.round(self.template['aspect_rate']* self.height).astype(int)


        self.panels = {int(index):value for index, value
                       in self.template.items() if index.isdigit()}
        # Rescale panels
        for _, panel in self.panels.items():
            panel['bbox']['x0'] = np.round(panel['bbox']['x0'] * self.width ).astype(int)
            panel['bbox']['x1'] = np.round(panel['bbox']['x1']* self.width).astype(int)
            panel['bbox']['y0'] = np.round(panel['bbox']['y0']* self.height).astype(int)
            panel['bbox']['y1'] = np.round(panel['bbox']['y1']* self.height).astype(int)
            # aspect ratio
            panel['ar'] = (panel['bbox']['x1'] - panel['bbox']['x0']) / (panel['bbox']['y1'] - panel['bbox']['y0'])
            # Initialize the associate image of the panel as None
            panel['associated_img'] = None
            panel['image_id'] = None

        # Create a template image, in which each pixels level represent a different class panel
        # Pixels with zeros values are background
        self.template_image = self._create_template_image()


    def _find_valid_template(self, query, dataset, templates, figure_type='inter-panel'):
        """
        This function select valid templates candidates for an image query.
        To consider a template valid it must have at least two panels
        with the aspect ratio of the input image, and all the others panels
        must be completed by a different image from the dataset.

        Parameters
        ----------

        query: dict or pandas dataframe
            {
                AR: <float>
                    Aspect Ratio of the image
                Class: <string>
                    Class of the image
            }

        dataset: pandas dataframe
            Columns:{
                    AR: <float>
                        Aspect Ratio
                    Class: <string>
                        Class of the image

            }
        templates: list(<string>)
            List of templates path

        figure_type: <string> in ['interPanel', 'intraPanel']
            Type of the compound figure

        """
        ar_image = query['AR']
        image_class = query['class']
        valid_templates = []
        if figure_type == 'inter-panel':
            panels_used_by_query = 2
        else:
            panels_used_by_query = 1

        # Check the templates that has more than one panel with the same AR as the input image
        for template_path in templates:
            with open(template_path,'r') as jf:
                template = json.load(jf)

            found = 0
            # check if we have data enough in the dataset to fill all panels
            image_allocated_to_figure = []
            for key, panel in template.items():
                if isinstance(panel,dict):

                    panel['bbox']['y1'] = y1  = np.round(panel['bbox']['y1'] * template['orig_height']).astype(int)
                    panel['bbox']['y0'] = y0  = np.round(panel['bbox']['y0'] * template['orig_height']).astype(int)
                    panel['bbox']['x1'] = x1  = np.round(panel['bbox']['x1'] * template['orig_width']).astype(int)
                    panel['bbox']['x0'] = x0  = np.round(panel['bbox']['x0'] * template['orig_width']).astype(int)

                    ar_panel = (x1-x0)/(y1-y0)
                    if (panel['class'] == image_class) and abs(ar_panel - ar_image) <= CompoundFigure.AR_TOLERENCE[image_class]:
                        found +=1

                if found >= panels_used_by_query:
                    break
            if found>=panels_used_by_query:
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

            # Reserve two panels to the forgery
            for index, row in candidate_panels.loc[(candidate_panels['class']==image_class) & ((ar_image-CompoundFigure.AR_TOLERENCE[image_class]) <= candidate_panels.AR) &
                                                (candidate_panels.AR <= (ar_image+CompoundFigure.AR_TOLERENCE[image_class]))].iloc[0:2].iterrows():
                candidate_panels.loc[index,'filled'] = 'query'

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


#-----------------------------------------------------------------------------------------------------#
#                                                                                                     #
#                           METHODS RELATED WITH IMAGE PANEL                                          #
#                                                                                                     #
#-----------------------------------------------------------------------------------------------------#

    def load_image(self, img_path, objs_mask_path=None):
        """
        Load an image from a path
        TODO: Lead with all kind of extensions problems that PIL Image
        could have
        """
        img = np.array(Image.open(img_path))
        if objs_mask_path:
            objs_map = np.array(Image.open(objs_mask_path).convert('L'))
        else:
            objs_map = np.zeros_like(img)
        return img, objs_map


    def _get_image_to_panel(self, panel, panel_class):
        """
        Choose a secundary role image that fits the panel_class, that does not participate in the forgery.
        This secundary image will be insert in the compound figure, but would participate from the forgery.

        """
        # Creates a list of candidate images that has the same class pointed in the template (panel_class) and fits
        # the panel
        #candidate_images = [ img_id for img_id, img in self._dataset[panel_class].items()
        #                        if np.absolute(img['img_ar']-panel['ar']) <= CompoundFigure.AR_TOLERENCE[panel_class]]
        cand_img = self._dataset.loc[(self._dataset['class']==panel_class) &
                                        ((self._dataset.AR - CompoundFigure.AR_TOLERENCE[panel_class]) <=panel['ar']) &
                                        (panel['ar'] <= self._dataset.AR + CompoundFigure.AR_TOLERENCE[panel_class])]
         # Choose a random image from the selection above
        cand_img = cand_img.sample(n=1).iloc[0]
        img_path = cand_img['dataPath']
        #image_id = random.choice(candidate_images)
        #candiate_image = self._dataset[panel_class][image_id]
        #img_path = candiate_image['img_path']

        # Delete the image from the dataset to avoid further contamination
        drop_content = self._dataset[self._dataset['dataPath'] == img_path].index
        self._dataset = self._dataset.drop(drop_content)
        # del self._dataset[panel_class][image_id]

        # Insert the panel image id with the used image
        panel['image_id'] = img_path

        image, _ = self.load_image(img_path)
        return image


    def _select_panel_id(self, panel_class, ar):
        """
        select a panel from self.panels, that is from the
        panel_class and have aspect ratio |p['ar']-ar| <= CompoundFigure.AR_TOLERENCE[panel_class]
        """
        candidate_panel = []
        for panel_id, panel in self.panels.items():
            if panel['class'] == panel_class and  \
               np.absolute( ar - panel['ar'] ) <= CompoundFigure.AR_TOLERENCE[panel_class] and \
               panel['associated_img'] is None:

                candidate_panel.append(panel_id)
        if len(candidate_panel) > 0:
            return random.choice(candidate_panel)

        else:
            assert False, f"None panel with the AR: {ar} and Class: {panel_class} remained"


#-----------------------------------------------------------------------------------------------------#
#                                                                                                     #
#                           METHODS RELATED WITH FIGURE ASSEMBLY                                      #
#                                                                                                     #
#-----------------------------------------------------------------------------------------------------#

    def _assembly_figure(self, forgery_panel_type):
        """
        Assembly all panels using the template and the images provided.
        A ground-truth figure will be generated as well. If the forgery is of type
        Intra-Panel panel, the  ground-truth will highlight just one panel; otherwise,
        it will highlight all panels involved.
        """

        # Create figure with a white background
        figure = np.ones((self.height,self.width,3))*255
        figure = figure.astype(np.uint8)

        # Generate and insert all graphs on the respective panels
        for _, panel in self.panels.items():
            if not panel['associated_img'] is None:
                panel_image = panel['associated_img']

            elif panel['class'] in self._dataset['class'].unique():
                panel_image = self._get_image_to_panel(panel, panel['class'])

            # panel is 'graph' or 'others'
            else:
                # Get shape of the panel
                figshape = (int(panel['bbox']['y1']) - int(panel['bbox']['y0']),
                            int(panel['bbox']['x1']) - int(panel['bbox']['x0']))
                # Randomly generate graph to insert in panel
                if np.random.randint(2):
                    panel_image = self.generate_fake_bar(figshape)
                else:
                    panel_image = self.generate_fake_plot(figshape)
                panel['image_id'] = 'RANDOM_GRAPH'

            figure = self.insert_panel_in_figure(figure, panel_image, panel['bbox'])

        self.figure = figure

        # Groundtruth
        self.figure_groundtruth_forgery = np.zeros((self.height,self.width,3))

        # If type of figure is not pristine
        if forgery_panel_type != "pristine-panel":

            forgery_panel = [p for _,p in self.panels.items() if p['class'] =='__FORGERY__'][0]
            self.figure_groundtruth_forgery = self.insert_panel_in_figure(self.figure_groundtruth_forgery,
                                                        self.forgery_image_gt,
                                                        forgery_panel['bbox'])
            self.figure_groundtruth_forgery = self.figure_groundtruth_forgery[:,:,0].astype(np.uint8)

            forgery_panel = [p for _,p in self.panels.items() if p['class'] =='__FORGERY__'][0]

            # Pristine ground-truth shows the objects involved in the forgery before
            # the manipulation was applied
            if forgery_panel_type == 'intra-panel':
                pristine_panel = forgery_panel
            else:
                pristine_panel = [p for _,p in self.panels.items() if p['class'] =='__PRISTINE__'][0]
            self.figure_groundtruth_pristine = np.zeros((self.height,self.width,3))
            self.figure_groundtruth_pristine = self.insert_panel_in_figure(self.figure_groundtruth_pristine,
                                                            self.pristine_gt,
                                                            pristine_panel['bbox'])
            self.figure_groundtruth_pristine = self.figure_groundtruth_pristine[:,:,0].astype(np.uint8)



#-----------------------------------------------------------------------------------------------------#
#                                                                                                     #
#                           METHODS RELATED WITH FORGERY EXECUTION                                    #
#                                                                                                     #
#-----------------------------------------------------------------------------------------------------#

    def _assert_panels_to_forgery(self, forgery_class, minum_number_of_panels=1):
            """
            Check if there is an available panel in the template layout that could fit the
            forgery_image
            """
            f_img, _ = self.load_image(self.forgery_image_path)
            # Image aspect ratio
            f_ar = f_img.shape[1] / f_img.shape[0]
            
            accepted_panels = 0
            for _, panel in self.panels.items():
                if panel['class'] == forgery_class and  \
                np.absolute( f_ar - panel['ar'] ) <= CompoundFigure.AR_TOLERENCE[forgery_class] and \
                    panel['associated_img'] is None:
                    accepted_panels+=1
                    if accepted_panels >= minum_number_of_panels:
                        return True

            assert False, "The forgery image could not fit the template"


    def _execute_inter_panel_forgery(self,forgery_info, forgery_class ):
            """
            This method creates a Inter-Panel forgery. It can uses the class SimpleFigureForgery
            to create a forgery image or call the method multi_image_forgery
            """
            function_name = forgery_info['function_name']
            # If the forgery function is intrisic multi panel, (copy a region from one panel to other)
            if function_name in CompoundFigure.MULTI_PANEL_FORGERY:
                    self._execute_intrinsic_multi_panel_forgery(forgery_info, forgery_class)

            # If the forgery function could be applied to only one panel, and
            # then copy pasted to another panel
            else:
                self._execute_intra_panel_forgery(forgery_info,forgery_class)

                # Select other panel with similar aspect rate and same forgery_class
                # to be the panel that we will copy from
                ar = self.pristine_image.shape[1] / self.pristine_image.shape[0]
                panel_id = self._select_panel_id(forgery_class, ar)
                self.panels[panel_id]['associated_img'] = self.pristine_image
                self.panels[panel_id]['class'] = '__PRISTINE__'
                # This is a case of panels copy-move so the pristine image path is the same as the forgery
                # path
                self.pristine_image_path = self.forgery_image_path
                self.panels[panel_id]['image_id'] = self.forgery_image_path

                # Fix Groundtruth, Since the same image will entirely be copy on two panels,
                # Both panels must represent this
                self.forgery_image_gt = np.ones_like(self.single_image_forgery.groundtruth)*255
                self.pristine_gt = np.ones_like(self.single_image_forgery.pristine_groundtruth)*255


    def _execute_intrinsic_multi_panel_forgery(self, forgery_info, forgery_class):
        """
        Apply forgery that only could be applied on a multi panel fashion
        """

        function_name = forgery_info['function_name']

        # Splincing function
        if 'splicing' in function_name:
            donor_img_path, donor_obj_path =  forgery_info['donor_img_path'], forgery_info['donor_obj_path']
            host_img_path, host_obj_path  = forgery_info['host_img_path'], forgery_info['host_obj_path']

            # Donor --> Prisitine Image
            self.pristine_image, _ = self.load_image(donor_img_path)
            # Apply forgery
            args = forgery_info['args']
            args['donor'], args['donor_map'] = self.load_image( donor_img_path, donor_obj_path)
            args['host'], args['host_map'] =  self.load_image( host_img_path, host_obj_path)
            forgery_call = f"{function_name}(**args)"
            self.forgery_image, self.pristine_gt, self.forgery_image_gt = eval(forgery_call)

            # Set images path to insert on tyhe panel id
            self.forgery_image_path = host_img_path
            self.pristine_image_path = donor_img_path


        # Forgery is a simple copy-move of a image
        else:

            # Check if there is a panel that could fit  the image aspect ratio
            self._assert_panels_to_forgery(forgery_class,2)

            # Apply forgery
            args = forgery_info['args']
            args['img'], _ = self.load_image( self.forgery_image_path)
            forgery_call = f"{function_name}(**args)"
            self.pristine_image, self.pristine_gt, self.forgery_image ,self.forgery_image_gt = eval(forgery_call)

            # Set images path to insert on tyhe panel id
            self.pristine_image_path = self.forgery_image_path

        # Associate forgery image with a panel
        ar = self.forgery_image.shape[1] / self.forgery_image.shape[0]
        panel_id = self._select_panel_id(forgery_class, ar)
        self.panels[panel_id]['associated_img'] = self.forgery_image
        self.panels[panel_id]['class'] = '__FORGERY__'
        self.panels[panel_id]['image_id'] = self.forgery_image_path
        # Associate pristine image with a panel
        panel_id = self._select_panel_id(forgery_class, ar)
        self.panels[panel_id]['associated_img'] = self.pristine_image
        self.panels[panel_id]['class'] = '__PRISTINE__'
        self.panels[panel_id]['image_id'] = self.pristine_image_path


    def _execute_intra_panel_forgery(self, forgery_info, forgery_class):

        # Check if there is a panel that could fit
        # the image aspect ratio
        self._assert_panels_to_forgery(forgery_class)

        # Execute forgery
        img, objs_map = self.forgery_image_path , self.forgery_obj_map_path
        self.pristine_image, _ = self.load_image(img)
        self.single_image_forgery = SimpleFigureForgery(img, forgery_info, objs_map)
        self.forgery_image = self.single_image_forgery.forgery_img
        self.forgery_image_gt = self.single_image_forgery.gt_after_forgery
        self.pristine_gt = self.single_image_forgery.gt_before_forgery

        # Associate forgery image with a panel
        ar = self.forgery_image.shape[1] / self.forgery_image.shape[0]
        panel_id = self._select_panel_id(forgery_class, ar)
        self.panels[panel_id]['associated_img'] = self.forgery_image
        self.panels[panel_id]['class'] = '__FORGERY__'
        self.panels[panel_id]['image_id'] = self.forgery_image_path

#-----------------------------------------------------------------------------------------------------#
#                                                                                                     #
#                               METHODS RELATED WITH METADATA                                         #
#                                                                                                     #
#-----------------------------------------------------------------------------------------------------#
    def _get_forgery_metadata(self, modality, figure_type):
        """
        Return the metada involved in the forgery
        """

        metadata = {}
        forgery_info_args = self.forgery_info.get('args').copy() if self.forgery_info.get('args') else dict()


        if self.forgery_info['function_name'] == 'splicing_forgery':
            del forgery_info_args['donor']
            del forgery_info_args['donor_map']
            del forgery_info_args['host']
            del forgery_info_args['host_map']
        else:
            if self.forgery_info['function_name'] != 'pristine':
                del forgery_info_args['img']
            if self.forgery_info.get('function_name') in CompoundFigure.OBJ_MAP_NEEDED:
                del forgery_info_args['objs_map']

        metadata['forgery_info'] = {
                'class': self.forgery_info['forgery_class'],
                'function_name': self.forgery_info['function_name'],
                'args': forgery_info_args,
                "figure_type": figure_type,
                'modality': modality,
        }
        # Forgery Panel annotation
        metadata["figure_annotations"] = {
                                "height": self.figure.shape[0],
                                "width": self.figure.shape[1],
                                "template":self.figure_template
                    }

        # Panels Annotation
        panels_annotations = copy.deepcopy(self.panels)
        for img_id, img_id_ann in panels_annotations.items():
            del panels_annotations[img_id]['associated_img']
            del panels_annotations[img_id]['ar']
            # Serializing object
            for key, val in panels_annotations[img_id]['bbox'].items():
                panels_annotations[img_id]['bbox'][key] = int(val)

        metadata['figure_annotations'].update(panels_annotations)

        return metadata





def insert_text( figure,
                text,
                x,y,
                font_name="forgery_lib/figure/fonts/times-new-roman-14.ttf",
                font_size=45,
                color_font='rgb(0,0,0)'):
    """
    Insert text in the figure in the location of bbox.
    This function is used only for caption, legend propouse

    Return
    Figure with the text inserted
    """


    fig = Image.fromarray(figure)
    draw = ImageDraw.Draw(fig)
    font = ImageFont.truetype(font_name, size=font_size)
    draw.text((x,y), text,
                    fill=color_font, font=font)

    return np.array(fig)

def get_color_font(figure,bbx):
    """
    Based on the histogram of the bbx, return a contrast color
    """

    figure_pil= Image.fromarray(figure)
    if not 'L' in figure_pil.getbands():
        figure_pil.convert("L")

    # Check histogram of bouding box
    img_array  = np.array(figure_pil)
    img_array = img_array[bbx]
    hist, bins = np.histogram(img_array)
    # If value of hist concentrate in a dark region
    if bins[np.argmax(hist)] < 100:
        return "rgb(255,255,255)"
    # If value of hist concentrate in bright region
    return "rgb(0,0,0)"
