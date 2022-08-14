
# Please cite the paper if you use this code:
# Y. M. Li and J. T. Zhou, ��Image Copy-Move Forgery Detection Using Hierarchical Feature Point Matching��, ASC, 2016.
#########################################################

import cv2
import numpy as np
import matplotlib.pyplot as plt

class Parameters(object):

	def __init__(self, img) -> None:

		h, w = img.shape[:2]
		self.max_sup_radius = 20
		self.min_sup_radius = 16
		self.max_hole_size = 300
		self.max_keypoints_octave = 7000
		self.r_ratio = 16
		self.check_orienation  = 1
		self.match_strategy = 1
		self.post_strategy = 1
		self.allowed_err = 1
		self.check_scale = 1
		self.check_dis = 0
		self.check_edges = 0
		self.check_transition= 1 #for coverage database
		self.illumination=0 #for coverage database
		self.max_scale = 2 #for coverage database
		self.fillholes =0 #for coverage database
		self.allow_error_inliers=1.5 #for coverage database
		self.half_dis_add = 50;
		self.scale_err_ratio = 5;
		self.match_ratio = 1/5; #the correct match at least 1/5 of total matches

		if w < 840 and h<840:
			self.check_orienation=1;
			self.match_strategy = 2;
			self.post_strategy=2;
			self.bi_direction = 1;
			self.sift_method =1;
			self.match_thr = 0.74;
			self.min_inliers = 4;
			self.min_clone_dis = max(40,25*(np.sqrt(w*h/(1000*1000))));
			self.scale_seg = [0];
			self.min_size = 250;
			self.min_total_inliers = 4;
			self.min_neighbors = min(3, int(np.round(3*(w*h/(1000*1000)))));
			if w<600 and h<600:        #for coverage database
				self.sift_method =1;
				self.check_transition=0;
				self.match_thr = 0.65;
				self.matching_rgb_dif = 70;
				self.min_size = 2000;
				self.min_inliers = 5;
				self.min_total_inliers =6;
				self.allow_error_inliers=1.5;
				self.illumination=15;
				self.fillholes =1;
				self.max_scale = 2.5;
		else:
			self.match_strategy = 2;
			self.bi_direction = 1;
			self.match_thr = 0.52;
			self.sift_method =1;
			self.min_inliers = 5;
			self.min_clone_dis = 25*(np.sqrt(w*h/(1000*1000)));
			self.scale_seg = [-1,0,1,2,3];
			self.min_size = min(500,400 + int(np.round(30*w*h/(1000*1000))));
			self.min_neighbors = min(3, int(np.round(3*np.sqrt(w*h/(1000*1000)))));
			self.min_total_inliers = 20;
			self.scale_err_ratio = 4;
			self.match_ratio = 1/5; #the correct match at least 1/5 of total matches
			self.half_dis_add = 50;
			self.min_dis = 25;
			self.check_dis = True;
			self.check_edges = True;
		if w>=1400 and h>=1400:
			self.match_strategy = 1;
			self.half_dis_add = 80;
			self.allowed_err = 8;
			self.bi_direction = 2;
			self.min_inliers = 6;  #5
			self.max_sup_radius = 48;
			self.min_sup_radius = 20;
			self.max_hole_size = 1500;
			self.max_keypoints_octave = 16000;
			self.min_total_inliers = 30;
			self.scale_err_ratio = 30;
			self.check_scale = 0;
			self.match_ratio =1/20;
			self.scale_seg = [-1,0,1];
			self.match_thr = 0.6;                  #for F660, set it as 0.5, can achieve much lower FPR
			self.min_neighbors = min(2, int(np.round(3*np.sqrt(w*h/(1000*1000)))));
			self.min_dis = 130;

		if self.bi_direction ==2:
			self.min_size = int(np.round(self.min_size*1.3));
		self.max_matches   = 50000;
		self.max_err_dis_xy = 100;

def process_image (imagefile, plotimg, large_resize=0):
	"""
	Return
		bool_temp, map, inliers1, inliers2
	"""

	# read image
	img_brg = cv2.imread(imagefile)
	img_rgb = img_brg[...,::-1]

	# set_parameters
	para = Parameters(img_rgb)
	para.plotimg = plotimg

	if large_resize == 1:    #for resize factor bigger than 2
		para.scale_seg = -1;


	############## do matching ###############
	if para.match_strategy == 1:
		para.p1, para.p2 =  feature_matching(imagefile, para);
	else:
		para.p1, para.p2 =  feature_matching2(imagefile, para);

	#############draw figure ###################
	#if plotimg==1 and para.p1.shape[1]>=1:
	#	draw_match(img_rgb, para.p1, para.p2)


	###########postprocessing%########
	bool_temp, map, inliers1, inliers2 = post_processing(img_rgb, para);
	"""
	if para.post_strategy==1:
		bool_temp, map, inliers1, inliers2 = post_processing(img_rgb, para);
	else:
	    bool_temp, map, inliers1, inliers2 =   post_processing2(img_rgb, para);
	"""
"""
def draw_match(img, p1, p2):
plt.figure()
plt.imshow(img);
h,w = img.shape[:2];
line_width = max(1, int(w/1024 * h/768));
line_width = min(line_width, 1);
for i in range(p1.shape[1]):)
    plt.plot( (p1[1,i], p2[1,i]), [p1(2,i)' p2(2,i)'], 'color', 'blue','LineWidth', line_width);
scatter(p1(1,:),p1(2,:),'r');
scatter(p2(1,:),p2(2,:),'r');

"""

