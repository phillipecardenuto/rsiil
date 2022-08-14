import cv2
import numpy as np
from scipy.spatial.distance import pdist, squareform
import cyvlfeat
from PIL import Image

def gray_cluster(pixel_vals, step, length ,min_pixels):
    """
    @Return clusters, pixel cls
    """
    # pixel_vals = im[np.round(idx)];
    # print("--"*10)
    cl = np.arange(0,256,step)[:, np.newaxis];
    cl = cl-length;
    cl = np.concatenate([cl, cl+step+length], axis=1);
    cl[0,0] =0; cl[-1,-1] = 255;
    num_cls = cl.shape[0];
    bin_s  = np.zeros((1, num_cls));
    clusters =  [];
    pixel_cls = cl;
    re_cls = [];
    for i in range(num_cls):
        cur_range = cl[i,:];
        cur_samples = np.where(  np.logical_and( cur_range[0] <= pixel_vals , pixel_vals <= cur_range[1]))[0]

        if (len(cur_samples) <min_pixels):
            re_cls.append(i);
            continue;
        clusters.append(cur_samples);

    pixel_cls = np.delete(pixel_cls, re_cls, axis=0)
    return clusters, pixel_cls





def feature_matching(imgfile, para):
    """
    @Return:
        p1, p2
    """

    match_thr = para.match_thr;
    #im_brg = cv2.imread(imgfile)
    im_rgb = Image.open(imgfile)
    # bgr to rgb
    #im_rgb = im_brg[...,::-1].astype(np.float64)
    im_rgb = np.array(im_rgb).astype(np.float64)
    h, w       = im_rgb.shape[:2];

    levels = 3;
    min_octave = 0;
    # Gray image
    if len(im_rgb.shape) > 2:
        im = cv2.cvtColor(im_rgb.astype(np.uint8), cv2.COLOR_RGB2GRAY).astype(np.float64);
    else:
        im = im_rgb
    average_kernel = np.ones((3,3)) / 9
    im_smooth =  cv2.filter2D(im, -1, average_kernel);
    if para.sift_method == 1:
        peak_Tresh = 0.1;
        locs, descs = cyvlfeat.sift.sift(
                                image = im.astype(np.float32),
                                n_levels = levels,
                                first_octave = min_octave,
                                n_octaves = 5,
                                peak_thresh = peak_Tresh,
                                edge_thresh = 12,
                                compute_descriptor = True)
                                #'Max_keypoints',para.max_keypoints_octave ) ;
        feature_len = locs.shape[0];
        #print('found %d features\n'%feature_len);
    elif para.sift_method == 2:
        raise NotImplementedError
        # sift = cv2.SIFT_create()
        # locs, descs = sift.detectAndCompute(im.astype(np.uint8), None)

    #### the min_sigma for each octave

    # locs -> [:, [Y,X,S,TH]]
    # where ``(Y, X)`` is the floating point center of the
    # keypoint, ``S`` is the scale and ``TH`` is the orientation (in radians).

    #locs = locs.transpose() # no need to transpose here
    locs[:, [1,0]] = locs[:, [0,1]] # matlab --> locs = locs(:,[2,1,3,4]); # change Y,X to X,Y
    # descs = descs.transpose().astype(np.float64) # no need to transpose here
    sigmak = 2**(1/levels);
    sigma_octave = sigmak**levels;
    sigma0 = 1.6;
    octaves_idx  = para.scale_seg;
    octaves_sigma = np.power(sigma0*sigma_octave, octaves_idx);
    ##############################################

    detected_sigmas = locs[:,2]; # get scales
    extend_octaves_sigma = np.concatenate([[-np.inf], octaves_sigma[1:], [np.inf]])
    locs_scales    =  [];
    descs_scales =  [];
    octaves_num  =  len(octaves_idx);
    for i in range(octaves_num):
        sigma_low = extend_octaves_sigma[octaves_num -i-1];
        sigma_high = extend_octaves_sigma[octaves_num -i];

        locs_scales.append(
            np.squeeze(
                locs[ np.where(np.logical_and(sigma_low<=detected_sigmas, detected_sigmas<sigma_high)),:]
            )
                );
        descs_scales.append(
            np.squeeze(
                descs[ np.where(np.logical_and(sigma_low<=detected_sigmas, detected_sigmas<sigma_high)),:]
            )
        )

    ##########do matching%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1=np.array([]); p2=np.array([]);
    num = 0;

    for oct_idx in range(len(octaves_idx)):
        cur_locs = locs_scales[oct_idx];
        cur_des =  descs_scales[oct_idx];
        if len(cur_locs) == 0:
            #print('no keypoints are detected in %dth octave\n'% (oct_idx+1));
            continue;
        if cur_locs.ndim == 1:
            cur_locs = np.expand_dims(cur_locs, axis=0)
        if cur_des.ndim == 1:
            cur_des = np.expand_dims(cur_des, axis=0)
        temp_locs = np.round(cur_locs[:,:4]).astype(np.int);
        temp_locs[:,2:4] = cur_locs[:,2:4];

    ##############gray level cluster%%%%%%%%%%%%%%%%%%%%%%
        cur_keys_n = temp_locs.shape[0]
        gray_cls = 1;
        clusters = [];
        if cur_keys_n >=5000:
            # this is the matrix addess just like in language C  M[x,y] = y * (M)width + x
            # key_indx = (temp_locs[:,1]-1)*w + temp_locs[:,0]; --> not work on python

            píxel_vals = im_smooth[temp_locs[:,1], temp_locs[:,0]]
            if cur_keys_n < 10000:
                clusters, _ = gray_cluster(píxel_vals,  40,10, 5);
            else:
                clusters, _ = gray_cluster(píxel_vals,  20,5, 5);
            gray_cls = len(clusters);

        else:
            clusters.append(np.arange(cur_keys_n));

        for gray_indx in range(gray_cls):
            loc1 = temp_locs[clusters[gray_indx], :];
            des1 = cur_des[clusters[gray_indx],:];
            #div_factor = (np.sqrt(np.sum(des1*des1, 1)))
            norm = np.linalg.norm(des1,axis=1)
            #des1 = np.divide(des1, div_factor[:, np.newaxis])
            des1 = des1 / norm[:,None] # Normalize each vector by its norm ( set des1/|des1| =1 )
            #MATLAB --> des1 = des1./repmat(sqrt(sum(des1.*des1,2)),1,size(des1,2));
            ### sift matching ###
            des2t = des1.transpose();
            match_indx = [];
            if des1.shape[0] > 1:
                for i in range( des1.shape[0]):
                    if i in match_indx:
                        continue;
                    dotprods = des1[i,:] @ des2t;
                    arc_cos = np.arccos(np.round(dotprods, decimals=7)) # round to avoid 1.0000000000000000000002
                    indx = np.argsort(arc_cos);
                    j=1;
                    if len(indx)<= j+1:
                        continue
                    if arc_cos[indx[j]] < match_thr* arc_cos[indx[j+1]]:
                        match_i = indx[j];
                        if np.linalg.norm( loc1[i,:]-loc1[match_i,:], 2) > para.min_clone_dis:
                            num=num+1;
                            loc1_i = loc1[i,:]
                            if len(loc1_i.shape) == 1:
                                loc1_i = np.expand_dims(loc1_i, axis=0)

                            if len(p1)==0:
                                p1 = np.squeeze(loc1_i)
                            else:
                                p1 = np.vstack([p1, np.squeeze(loc1_i)])

                            loc1_matchi = loc1[match_i,:]
                            if len(loc1_matchi.shape) == 1:
                                loc1_matchi = np.expand_dims(loc1_matchi, axis=0)
                            if len(p2) == 0:
                                p2 = np.squeeze(loc1_matchi)
                            else:
                                p2 = np.vstack([p2, np.squeeze(loc1_matchi)])
                            match_indx.append( (i,match_i));
        if len(p1):
            if p1.ndim == 1:
                p1 = np.expand_dims(p1, axis=0)
                p2 = np.expand_dims(p2, axis=0)
            p = np.round( np.hstack( [p1[:,:2], p2[:, :2]]));
            p_temp, indx = np.unique(p, return_index=True, axis=0);
            p1 = np.hstack([p_temp[:,:2], p1[indx,2:]]);
            p2 = np.hstack([p_temp[:,2:], p2[indx,2:]]);
            num=p1.shape[0]
    #print('Found %d matches.\n'%num);

    ###################remove isolated matching%%%%%%%%%%%%%%%%%
    num = 0;
    if len(p1):
        half_dis = (para.min_clone_dis/2) + para.half_dis_add;
        dis_map = pdist( np.vstack([p1[:,0:2], p2[:,0:2]]));
        dis_map_q = squareform(dis_map <= half_dis);
        neighbors = np.sum(dis_map_q, axis=0);
        neighbors_p1 = neighbors[:len(p1)];
        neighbors_p2 = neighbors[len(p1):];
        indx = np.logical_or(neighbors_p1>=para.min_neighbors, neighbors_p2 >=para.min_neighbors);
        num = np.sum(indx, axis=0);
        p1 = np.vstack([ np.transpose(p1[indx,0:2]), np.ones((1, num)) , np.transpose(p1[indx,2:])]);
        p2 = np.vstack([ np.transpose(p2[indx,0:2]), np.ones((1, num)) , np.transpose(p2[indx,2:])]);
    #print('after remove the isoloated matching, found %d matches\n'% num);

    # First line of p1 contains all points x coordinate
    # Second line of p1 contains all points y coordinate
    # Third contains all ones
    # Fourth contains the scales
    # Fifth contains the thresholds
    return p1, p2
