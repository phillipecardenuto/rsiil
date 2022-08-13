# 'require' this file in your command generator

# extend or modify this file to create your own configurations

# nul: directly combines the images, using the preset coordinates
# rot: applies rotation in steps of 10 degrees up to 180 degrees, adds some noise
# rot_noise: applies rotation in steps of 10 degrees up to 180 degrees, adds
#            different levels of global noise
# too_much_dont_try: iterates over all available parameters, mainly for
#                    demonstration purposes. Note that every input image
#                    creates 7600 benchmark images, which is typically a bit
#                    too much for evaluation.

our %splice_config = (
	nul => {
		gt_type => 'cmfd', # alternatively: gt_type => 'snippet'
	},

	rot => { # resizes images w/ nearest neighbor interpolation
		rot_min => 2, # 5 cases
		rot_max => 10,
		rot_step => 2,
		gt_type => 'cmfd',
	},
	rot_noise => {
		# 1*19*6 = 114 cases per input image; already pretty much :/
		lnoise_min => 0, # 1 case
		lnoise_max => 0,
		lnoise_step => 0,
		rot_min => 0, # 19 cases
		rot_max => 180,
		rot_step => 10,
		gnoise_min => 0, # 6 cases
		gnoise_max => 50,
		gnoise_step => 10,
		gt_type => 'cmfd',
	},
	scale_down => { # resizes images w/ nearest neighbor interpolation
		cont_scale_min => 100, # 5 cases
		cont_scale_max => 1000,
		cont_scale_step => 200,
		gt_type => 'cmfd',
	},
	cmb_easy1 => {
		scale_min => 1010, # 1 case
		scale_max => 1010,
		scale_step => 100,
		rot_min => 2, # 1 case
		rot_max => 2,
		rot_step => 1,
		gjpeg_min => 80, # 1 case
		gjpeg_max => 80,
		gjpeg_step => 10,
		gt_type => 'cmfd', # alternatively: gt_type => 'snippet'
	},
	orig_sd => {
		cont_scale_min => 500, # 5 cases
		cont_scale_max => 500,
		cont_scale_step => 100,
		gt_type => 'none', # alternatively: gt_type => 'snippet'
	},
	orig_jpeg_sd => {
		cont_scale_min => 500, # 5 cases
		cont_scale_max => 500,
		cont_scale_step => 100,
		gjpeg_min => 20, # 5 cases
		gjpeg_max => 100,
		gjpeg_step => 10,
		gt_type => 'none', # alternatively: gt_type => 'snippet'
	},
	orig_jpeg => { # resizes images w/ nearest neighbor interpolation
		gjpeg_min => 20, # 5 cases
		gjpeg_max => 100,
		gjpeg_step => 10,
		gt_type => 'none', # alternatively: gt_type => 'snippet'
	},
	jpeg => { # resizes images w/ nearest neighbor interpolation
		# -> per input image, 9 output images are created
		gjpeg_min => 20, # 5 cases
		gjpeg_max => 100,
		gjpeg_step => 10,
		gt_type => 'cmfd', # alternatively: gt_type => 'snippet'
	},
	lnoise => { # resizes images w/ nearest neighbor interpolation
		# -> per input image, 9 output images are created
		lnoise_min => 20, # 5 cases
		lnoise_max => 100,
		lnoise_step => 20,
		gt_type => 'cmfd', # alternatively: gt_type => 'snippet'
	},
	scale => { # resizes images w/ nearest neighbor interpolation
		scale_min => 910, # 10 cases
		scale_max => 1090,
		scale_step => 20,
		gt_type => 'cmfd',
	},
	nul_sd => {
		cont_scale_min => 500, # 5 cases
		cont_scale_max => 500,
		cont_scale_step => 100,
		gt_type => 'cmfd', # alternatively: gt_type => 'snippet'
	},
	jpeg_sd => { # resizes images w/ nearest neighbor interpolation
		cont_scale_min => 500, # 5 cases
		cont_scale_max => 500,
		cont_scale_step => 100,
		gt_type => 'cmfd',
		# -> per input image, 9 output images are created
		gjpeg_min => 20, # 5 cases
		gjpeg_max => 100,
		gjpeg_step => 10,
		gt_type => 'cmfd', # alternatively: gt_type => 'snippet'
	},
	rot_sd => { # resizes images w/ nearest neighbor interpolation
		cont_scale_min => 500, # 1 cases
		cont_scale_max => 500,
		cont_scale_step => 100,
		rot_min => 2, # 5 cases
		rot_max => 10,
		rot_step => 2,
		gt_type => 'cmfd',
	},
	scale_sd => { # resizes images w/ nearest neighbor interpolation
		scale_min => 910, # 10 cases
		scale_max => 1090,
		scale_step => 20,
		cont_scale_min => 500, # 1 cases
		cont_scale_max => 500,
		cont_scale_step => 100,
		gt_type => 'cmfd',
	},
	resampling => {
		scale_min => 500, # 5 cases
		scale_max => 1500,
		scale_step => 100,
		gt_type => 'snippet', # alternatively: gt_type => 'cmfd' or 'none'
	},
	resampling2 => {
		scale_min => 800, # 5 cases
		scale_max => 1500,
		scale_step => 600,
		gjpeg_min => 80,
		gjpeg_max => 100,
		gjpeg_step => 10,
		gt_type => 'snippet', # alternatively: gt_type => 'cmfd' or 'none'
	},
	too_much_dont_try => {
		# demo of the possible parameters. Never execute it, it's a waste of
		# hard disk space.
		# -> per input image, 4*19*5*5*4=7600 output images are created => the
		# whole hard drive is full of trash; plus the evaluation takes forever
		lnoise_min => 5, # 4 cases
		lnoise_max => 20,
		lnoise_step => 5,
		disc_scale_min => 850, # 6 cases
		disc_scale_max => 1100,
		disc_scale_step => 50,
		cont_scale_min => 850, # 6 cases
		cont_scale_max => 1100,
		cont_scale_step => 50,
		rot_min => 0, # 19 cases
		rot_max => 180,
		rot_step => 10,
		scale_min => 800, # 5 cases
		scale_max => 1200,
		scale_step => 100,
		ljpeg_min => 55, # 5 cases
		ljpeg_max => 95,
		ljpeg_step => 10,
		gjpeg_min => 55, # 5 cases
		gjpeg_max => 95,
		gjpeg_step => 10,
		gnoise_min => 5, # 4 cases
		gnoise_max => 20,
		gnoise_step => 5,
		gt_type => 'cmfd', # alternatively: gt_type => 'snippet'
	},
);


