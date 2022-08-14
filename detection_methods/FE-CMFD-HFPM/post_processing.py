import cv2
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skimage import measure

def generated_object_match_map(map, p1, p2):
    if len(p1) == 0:
        return map

    # Remove all keypoints that both p1 and p2 do not fall in a region
    remove_points = []
    for p_index in range(p1.shape[1]):
        point1 = (int(p1[1,p_index]), int(p1[0,p_index]))
        point2 = (int(p2[1,p_index]), int(p2[0,p_index]))
        if map[point1] ==0 or map[point2] == 0:
            remove_points.append(p_index)

    p1_consistent = np.delete(p1, remove_points, axis=1)
    p2_consistent = np.delete(p2, remove_points, axis=1)

    # Label connected components
    label_map, num_conn_comp = measure.label(map.astype(int), background=0, return_num=True)

    # associate each keypoint with a connected component
    cc_p1 = []
    cc_p2 = []
    connected_components_found = []
    for p_index in range(p1_consistent.shape[1]):
        point1 = p1_consistent[:2, p_index]
        point2 = p2_consistent[:2, p_index]
        # fix x and y position
        point1 = point1[::-1].astype(int)
        point2 = point2[::-1].astype(int)

        found1 = False
        found2 = False
        insert1 = False
        insert2 = False

        for cc in range(1, num_conn_comp+1):
            cc_map = np.argwhere(label_map == cc)

            if point1 in cc_map:
                if not label_map[tuple(point1)] in connected_components_found:
                    cc_p1.append((point1,label_map[tuple(point1)], p_index+1)) # point coord, connected component coord, point label
                    connected_components_found.append(label_map[tuple(point1)])
                    insert1 = True
                found1 = True

            if point2 in cc_map:
                if not label_map[tuple(point2)] in connected_components_found:
                    cc_p2.append((point2, label_map[tuple(point2)], p_index+1)) # point coord, connected component coord, point label
                    connected_components_found.append(label_map[tuple(point2)])
                    insert2 = True
                found2 = True
            if found1 and  found2:
                break
        if insert1 ^ insert2: # xor
            if insert2:
                cc_p1.append((point1, label_map[tuple(point1)], p_index+1))
            else:
                cc_p2.append((point2, label_map[tuple(point2)], p_index+1))

    # Label final map
    final_label_map = - label_map.copy()
    index = 1
    while len(cc_p1) >0:
        region1 = cc_p1.pop()
        region2 = cc_p2.pop()
        conquired = False
        # If regions already conquired, assign the same label to the region2 (matching of region1)
        if (final_label_map[label_map == region1[1]] >0).any():
            final_label_map[label_map==region2[1]] = final_label_map[label_map == region1[1]][0]
            conquired = True

        # If regions already conquired, assign the same label to the region1 (matching of region2)
        if (final_label_map[label_map == region2[1]] >0).any():
            final_label_map[label_map==region1[1]] = final_label_map[label_map == region2[1]][0]
            conquired = True

        if not conquired:
            final_label_map[label_map==region1[1]] = index
            final_label_map[label_map==region2[1]] = index
            index += 1

    final_label_map[final_label_map<=0] =0
    return final_label_map

def hnormalise(x):
    """
    @ Return nx
    """

    rows,npts = x.shape;
    nx = x.copy();

    # Find the indices of the points that are not at infinity
    finiteind = np.where(abs(x[rows-1,:]) > np.finfo(float).eps)[0];

    #if len(finiteind) != npts:
        #%warning('Some points are at infinity');
        #print('Some points are at infinity');

    # Normalize points not at infinity
    for r in range(rows-1):
        nx[r,finiteind] = x[r,finiteind] / x[rows-1,finiteind];

    nx[rows-1,finiteind] = 1;
    return nx

def check_orientation(cos_t, sin_t, o1, o2, t):
    """
    @Return inliers
    """
    cos_t = np.round(cos_t, decimals=5) # Avoid 1.0000002
    theta = np.arccos(cos_t);
    if (cos_t<=0 and sin_t<=0) or (cos_t>0 and sin_t<0):
        theta = 2*np.pi - theta;

    est_o = np.round((theta/np.pi)*180, decimals=1);

    p1_o = np.round((o1 / np.pi)*180 , decimals=1)
    p2_o = np.round((o2 / np.pi)*180, decimals=1)
    p1_o[p1_o>=360] -= 360
    p2_o[p2_o>=360] -= 360

    p1_o[p1_o<0] = p1_o[p1_o<0] + 360;
    p2_o[p2_o<0] = p2_o[p2_o<0] + 360;

    dif_o = p2_o-p1_o;
    dif_o[dif_o<0] = dif_o[dif_o<0] + 360;
    dif_o2 = 360-dif_o;

    inliers1= np.abs(dif_o - est_o)<t;
    inliers2= np.abs(dif_o2 - est_o)<t;
    inliers = np.logical_or(inliers1, inliers2);
    return inliers

def noise_estimation(img, p1, p2):
    """
    @Return noise_err
    """
    h,w = img.shape[:2];

    s = 16;
    indx = np.logical_and(
            np.logical_and(
             np.logical_and(p1[0,:]>s, p1[1,:]>s),
             np.logical_and(p1[1,:]<w-s,p1[1,:]<h-s)
            ),
        np.logical_and(
            np.logical_and(
            p2[0,:]>s, p2[1,:]>s),
            np.logical_and(p2[1,:]<w-s, p2[1,:]<h-s)
        )
    )

    p1_ne = p1[:,indx].astype(np.int); p2_ne=p2[:,indx].astype(np.int);
    N = p1.shape[1];
    record_ratio = np.ones(N);
    for i in range(N):
        try:
            patch_1 = img[ p1_ne[0,i]-s:p1_ne[0,i]+s, p1_ne[1,i]-s:p1_ne[1,i]+s];
            patch_2 = img[p2_ne[0,i]-s:p2_ne[0,i]+s, p2_ne[1,i]-s:p2_ne[1,i]+s];
            [grad_tr_x1, grad_tr_y1] = np.gradient(patch_1.astype(np.float64), edge_order=2);
            mag_tr1 = np.mean(np.sum(np.sum(np.sqrt(grad_tr_x1**2+grad_tr_y1**2))));
            [grad_tr_x2, grad_tr_y2] = np.gradient(patch_2.astype(np.float64), edge_order=2);
            mag_tr2 = np.mean(np.sum(np.sum(np.sqrt(grad_tr_x2**2+grad_tr_y2**2))));
            record_ratio[i] = max(mag_tr1, mag_tr2)/max(1,min(mag_tr1, mag_tr2));
        except:
            continue

    if len(record_ratio) == 0:
        return 100
    record_ratio = np.sort(record_ratio);
    noise_ratio = np.mean(record_ratio[1:-1]);
    noise_err = 0;
    if noise_ratio > 10000:
        noise_err = 100;
    elif noise_ratio > 5000:
        noise_err = 60;
    elif noise_ratio > 1000:
        noise_err = 30;

    return noise_err

def check_distance(p1, p2, thr):
    """
    @Return
        check_pass
    """
    check_pass = True;
    abs_dis = abs(p1-p2);
    dis = np.sqrt((abs_dis[0,:]**2 )+(abs_dis[1,:]**2));
    if np.mean(dis) < thr:
        check_pass = False;
    return check_pass

def check_vh_edges(p1, p2, thr):
    """
    @Return
        check_pass
    """
    check_pass = True;
    N = p1.shape[1];
    temp1 = np.min([np.var(p1[0,:]), np.var(p1[1,:]), np.var(p2[0,:]), np.var(p2[1,:])]); #an vertical or horizontal edge
    abs_p1p2 = np.abs(p1-p2);
    temp2 = np.max([np.var(abs_p1p2[0,:]), np.var(abs_p1p2[1,:])]);

    if temp1<thr and temp2<1:
        check_pass = False;

    return check_pass

def get_inliers(H, p1, p2, t):
    """
    @ Return inliers
    """
    Hp1 = H @ p1;
    """
    Python could not successfully apply the inverse
    invHp2 = np.linalg.inv(Hp1) @ p2
    therefore, we are only cheking the H @ p1 towards p2
    """
    invHp2 = np.linalg.lstsq(H ,p2, rcond=1e-10 )[0];

    p1h     = hnormalise(p1.copy());
    p2h     = hnormalise(p2.copy());
    Hp1    = hnormalise(Hp1);
    invHp2 = hnormalise(invHp2);
    d2 = np.sum((p1h-invHp2)**2, axis=0)  + np.sum((p2h-Hp1)**2, axis=0);
    d2 = np.sum((p2-Hp1)**2, axis=0);
    inliers = abs(d2) < t;
    return inliers



def  post_processing(img_rgb, para):
    # Return [bool_temp, map, inliers1, inliers2]

    if len(img_rgb.shape) > 2:
        img = cv2.cvtColor(img_rgb.astype(np.uint8), cv2.COLOR_RGB2GRAY)
    else:
        img = img_rgb
    para.r_ratio = 16;
    r_ratio = para.r_ratio;
    p1 = para.p1.copy(); p2 = para.p2.copy();
    imclose_para = 5;
    max_sup_radius= para.max_sup_radius;

    h,w = img_rgb.shape[:2];

    # Return variables
    inliers1 = [];
    inliers2 = [];
    bool_temp=0;
    map = np.zeros((h,w));

    if len(p1) ==0:
        return bool_temp, map, inliers1, inliers2

    if  p1.shape[1] >= 5:

        average_kernel = np.ones((3,3)) / 9
        if len(img_rgb.shape) > 2:
            img_r =  cv2.filter2D(img_rgb[:,:,0].astype(np.float64), -1, average_kernel);
            img_g =  cv2.filter2D(img_rgb[:,:,1].astype(np.float64), -1, average_kernel);
            img_b =  cv2.filter2D(img_rgb[:,:,2].astype(np.float64), -1, average_kernel);
        else:
            img_gray = cv2.filter2D(img.astype(np.float64), -1, average_kernel);

        map_1 = np.zeros((h,w));
        map_2  = np.zeros((h,w));
        in_image_idx = np.logical_and(
                            np.logical_and(p1[0,:]<w, p1[1,:]<h),
                            np.logical_and(p2[0,:]<w, p2[1,:]<h)
        )
        p1 = p1[:,in_image_idx];
        p1[:, [1,0]] = p1[:, [0,1]]
        p2 = p2[:,in_image_idx];
        p2[:, [1,0]] = p2[:, [0,1]]
        p = np.hstack([p1[:2,:], p2[:2, :]]).T
        n_matches = p1.shape[1];
        if n_matches < 4:
            return bool_temp, map, inliers1, inliers2

        n_ranc = max(15, np.round(n_matches/30));
        n_ranc = min(n_ranc, max(10,np.round(n_matches/30)),40, n_matches);

        half_dis = (para.min_clone_dis/2)+55;
        dis_map = pdist(p)

        #@ Notes from the translator
        # Paper claims that half_dis = 100
        # But we are following the implemented solution


        dis_map_q = squareform(dis_map  <= half_dis);
        dis_map_q = np.logical_or(dis_map_q, np.eye(dis_map_q.shape[0]));


        indx = np.random.permutation(n_matches)
        n_ranc = int(n_ranc)
        indx = indx[:n_ranc];
        min_size = para.min_size;
        min_inliers = para.min_inliers;
        for k in range(n_ranc):
            cur_p = indx[k];
            neighbor_cur_p = dis_map_q[:,cur_p]
            cur_indx1 = neighbor_cur_p[:n_matches]
            cur_indx2 = neighbor_cur_p[n_matches:]


            or_cur_idx = np.logical_or(cur_indx1, cur_indx2)
            cur_p1 = p1[:, or_cur_idx];
            cur_p2 = p2[:, or_cur_idx];

            # The original implementation allows
            # the homography transformation to not include a
            # pair of matching keypoints, we did not follow this
            # strategy, since it appears differently in the paper

            # cur_p1 = p1[:, cur_indx1];
            # cur_p2 = p2[:, cur_indx2];

            # cur_indx_changed = np.logical_xor(
                            #  np.logical_or(cur_indx1, cur_indx2),
                            #   cur_indx1);

            # if np.sum(cur_indx_changed)>0:
                # cur_p1_add = p2[:, cur_indx_changed];
                # cur_p2_add = p1[:, cur_indx_changed];
                # cur_p1 = np.hstack([cur_p1, cur_p1_add]);
                # cur_p2 = np.hstack([cur_p2, cur_p2_add]);

            #######do ransac%%%%%%%%%%%%%%%%%
            if cur_p1.shape[1] < 4:
                #print('find matches less 4 in a local region, continue!\n ');
                continue;

            t  = 0.05;
            #[H, inliers1] = ransacfithomography2(cur_p1(1:3,:), cur_p2(1:3,:), t);
            min_local_size = min(cur_p1.shape[1], cur_p2.shape[1])
            # We are using homogeneous coordinates
            H, mask_inliers1 = cv2.findHomography(srcPoints=cur_p1[:3,:min_local_size].T,
                                             dstPoints=cur_p2[:3,:min_local_size].T,
                                             method=cv2.RANSAC,
                                             ransacReprojThreshold = 5,
                                             confidence=0.5)

            inliers1 = np.nonzero(mask_inliers1)[0] # get all inliers
            if (not np.array(H).any()) or len(inliers1)<4: #para.min_inliers
                #print('the number of inliers is small, continue!\n ');
                continue;

            ###### analyze the tranformation matrix (scale)%%%%%%%%%%%%%%%%
            allowed_err = para.allowed_err;
            try:
                [U,S,V] = np.linalg.svd(H[:2,:2]);
            except:
                return [bool_temp, map, inliers1, inliers2]

            if para.check_scale and (S[0]>8 or S[0]<0.2 or S[1]>8 or S[1]<0.2):
                #%print('estimated H ��scale�� is not accurate, give up\n');
                continue;
            cur_p1 = cur_p1[:,inliers1]
            cur_p2 = cur_p2[:,inliers1]

            rotation_xy = (U @ V).T;
            scale_err = max(abs(S[0]-1),abs(S[1])-1)*para.scale_err_ratio;
            coss = np.round(rotation_xy[0,0], decimals=7) # avoid coss of 1.0000000000002
            r_xy = (np.arccos(coss) /np.pi)*180;
            rotation_err = max(0,r_xy-5)*0.2;
            allowed_err = max(allowed_err,min(30,allowed_err + np.floor(scale_err+ rotation_err))); #50
            if not ((r_xy<30 or r_xy>70) and (r_xy<105 or r_xy>330)):
                min_inliers = max(4, para.min_inliers-1);
            if  para.check_orienation:
                indx_o  = check_orientation(
                                rotation_xy[0,0],
                                rotation_xy[1,0],
                                cur_p1[4,:],
                                cur_p2[4,:],
                                15);
                if np.sum(indx_o)<len(indx_o)*0.8:
                    #print('check dominant orientation fail, continue!\n ');
                    continue;

                # if sum(indx_o)<length(indx_o) && sum(indx_o)>4
                    # cur_p1 = cur_p1(:,indx_o); cur_p2 = cur_p2(:,indx_o);
                    # [H, inliers1] = ransacfithomography2(cur_p1(1:3,:), cur_p2(1:3,:), t);
                # end

            #do noise estimation only for those light attacks�� only
            #work for smooth images
            if (S[0]<1.2 and S[0])>0.8 and S[1]<1.2 \
                and S[1]>0.8 and (r_xy<30 or r_xy>330):
                noise_err = noise_estimation(img, cur_p1[0:2,:],cur_p2[0:2,:]);
                allowed_err = allowed_err + noise_err;
                if allowed_err>60:
                    min_inliers=4;

            ####### get all the inliers%%%%%%%%%%%%%%%%%%%%%%
            inliers =  get_inliers(H, p1[0:3,:], p2[0:3,:],allowed_err);
            outliers = np.where(inliers ==0)[0];

            cur_x1 = p1[:, inliers];
            cur_x2 = p2[:, inliers];

            if len(outliers) > 0:
                p1_1 = p2[:, outliers].copy();
                p2_1 = p1[:, outliers].copy();
                inliers_1 =  get_inliers(H, p1_1[0:3,:], p2_1[0:3,:], allowed_err);
                if np.sum(inliers_1)>0:
                    cur_x1 = np.hstack([cur_x1, p1_1[:, inliers_1]]);
                    cur_x2 = np.hstack([cur_x2, p2_1[:,inliers_1]]);

            if  para.check_orienation and len(cur_x1)>0:
                indx_o  = check_orientation(rotation_xy[0,0],
                                            rotation_xy[1,0],
                                            cur_x1[4,:],
                                            cur_x2[4,:],
                                            20);
                cur_x1 = cur_x1[:, indx_o];
                cur_x2 = cur_x2[:,indx_o];

                # if sum(~indx_o)>0
                    # fprintf('kick out some outliers through dominant orientation!\n ');
            if (cur_x1.shape[1]< max(min_inliers,np.floor(n_matches*para.match_ratio))) \
                    and cur_x1.shape[1]  < para.min_total_inliers:
                #print('detected random match pairs, give up!\n');
                continue

            if cur_x1.shape[1]<17 and para.check_dis and \
                 (not check_distance(cur_p1[0:2,:], cur_p2[0:2,:], para.min_dis)):
                #fprintf('check smaller distance fail, continue!\n ');
                continue;

            if cur_x1.shape[1]<7 and para.check_edges and \
                 allowed_err==para.allowed_err and  \
                 (not check_vh_edges(cur_x1[0:2,:],cur_x2[0:2,:], 25)):

                #%fprintf('suspect vertical or horizontal edge, give up!\n');
                continue;

            if S[0] > 1.1 or S[0] < 0.9 or S[1] > 1.1 or \
                S[1]<0.9:
                imclose_para = 7;
                max_sup_radius = para.max_sup_radius+12;
                if S[0]>2.8 or S[1]>2.8:
                    imclose_para = 15;
                    max_sup_radius = para.max_sup_radius+50;
                    r_ratio= r_ratio*2;

            for i in range(para.bi_direction):     #bi-direction transform



                H_cur=H.copy();
                if i == 1:
                    H_cur = np.linalg.inv(H);
                    temp = cur_x1.copy();
                    cur_x1 = cur_x2.copy();
                    cur_x2 = temp;
                #the pixels difference%%%%%%%%%%%%%%%%%%%%%%
                map_ori = np.zeros((h,w));
                cur_radius = cur_x1[3,:]*r_ratio;
                cur_radius[cur_radius>para.max_sup_radius] = max_sup_radius;
                cur_radius[cur_radius<para.min_sup_radius] = para.min_sup_radius;
                circles = np.transpose(np.vstack([cur_x1[0,:],
                                          cur_x1[1,:],
                                          cur_radius])).astype(np.int32)

                for circle in circles:
                    x,y, radius = circle
                    map_ori = cv2.circle(map_ori, (x,y), radius, 255, -1)
                map_ori[map_ori>0] = 1

                max_error = 16;#%;
                sus_region_1 =  np.array(np.where(map_ori ==1))
                r_indx = sus_region_1[0]  #np.floor((sus_region_1-1)/w)+1;
                c_indx = sus_region_1[1] #  sus_region_1 - (r_indx-1)*w;
                copy_region = H_cur @ np.vstack([r_indx, c_indx,np.ones_like(r_indx)])# np.ones(len(r_indx))])
                copy_region = np.round(copy_region[:2,:])

                sus_region_2 = copy_region.copy()
                sus_region_1 = sus_region_1[:,np.logical_and(
                                    np.logical_and(0<=sus_region_2[0], sus_region_2[0]<h),
                                    np.logical_and( 0<=sus_region_2[1], sus_region_2[1]<w) )].astype(np.uint)
                sus_region_2 = sus_region_2[:,np.logical_and(
                                    np.logical_and(0<=sus_region_2[0], sus_region_2[0]<h),
                                    np.logical_and( 0<=sus_region_2[1], sus_region_2[1]<w) )].astype(np.uint)

                if len(img_rgb.shape) > 2:
                    dif_r = abs(img_r[tuple(sus_region_1)] - img_r[tuple(sus_region_2)]);
                    dif_g = abs(img_g[tuple(sus_region_1)] - img_g[tuple(sus_region_2)]);
                    dif_b = abs(img_b[tuple(sus_region_1)] - img_b[tuple(sus_region_2)]);

                    refind_indx = np.logical_and(np.logical_and(dif_r< max_error, dif_g<max_error)
                                             ,dif_b<max_error);

                else:
                    dif = abs(img_gray[tuple(sus_region_1)] - img_gray[tuple(sus_region_2)]);
                    refind_indx = dif <max_error

                map_ori[tuple(sus_region_1[:, np.logical_not(refind_indx)])] = 0;


                map_copymove = np.zeros((h,w));
                map_copymove[tuple(sus_region_2[:, refind_indx])]=1;
                if i==1:
                    temp = map_ori;
                    map_ori = map_copymove;
                    map_copymove = temp;
                map_1 = np.logical_or(map_1, map_ori);
                map_2 = np.logical_or(map_2, map_copymove);


        if np.sum(map_1)>0 and np.sum(map_2)>0 \
            and (S[0]<4 or S[1]<4):
            map_1 =  bwareaopen(map_1, min_size);
            map_2 =  bwareaopen(map_2, min_size);

        if np.sum(map_1)>0 and np.sum(map_2)>0:
            map = np.logical_or(map_1, map_2);
            kernel_disk = get_kernel_disk(imclose_para)
            map = map.astype("uint8")
            map= cv2.morphologyEx(map, cv2.MORPH_CLOSE, kernel_disk)
            map = fill_small_holes(map, para.max_hole_size);
            bool_temp=1;
        if np.sum(map_1)>0 and np.sum(map_2)>0 and (S[0] >=4 or S[1]>=4):
            map =  bwareaopen(map, min_size);

    #if(bool_temp):
        #print('Tampering detected!\n\n');
    #else:
        #print('Image not tampered.\n\n');
    return [bool_temp, map, inliers1, inliers2]

def fill_small_holes(map, N):
    """
    @Return map_new
    """
    filled_map = bwfill(map);
    filled_map[filled_map>0] = 1
    holes = np.logical_and(filled_map, np.logical_not(map));
    bigholes = bwareaopen(holes, N);
    smallholes = np.logical_and(holes, np.logical_not(bigholes));
    map_new = np.logical_or(map, smallholes);
    return map_new

def bwfill(im_in):
    """
    Source: https://learnopencv.com/filling-holes-in-an-image-using-opencv-python-c/
    """

    th, im_th = cv2.threshold(im_in, 1, 255, cv2.THRESH_BINARY_INV);

    # Copy the thresholded image.
    im_floodfill = im_th.copy()

    # Mask used to flood filling.
    # Notice the size needs to be 2 pixels than the image.
    h, w = im_th.shape[:2]
    mask = np.zeros((h+2, w+2), np.uint8)

    # Floodfill from point (0, 0)
    cv2.floodFill(im_floodfill, mask, (0,0), 255);

    # Invert floodfilled image
    im_floodfill_inv = cv2.bitwise_not(im_floodfill)

    # Combine the two images to get the foreground.
    im_out = im_th | im_floodfill_inv
    return im_out


def bwareaopen(img_map, min_size):
    """
    Source: https://stackoverflow.com/a/42812226
    """
    img_map = img_map.astype("uint8")
    #find all your connected components
    nb_components, output, stats, _ = cv2.connectedComponentsWithStats(img_map, connectivity=4)
    #connectedComponentswithStats yields every separated component with information on each of them, such as size
    #the following part is just taking out the background which is also considered a component, but most of the time we don't want that.
    sizes = stats[1:, -1]; nb_components = nb_components - 1

    result = np.zeros((output.shape))
    #for every component in the image, you keep it only if it's above min_size
    for i in range(0, nb_components):
        if sizes[i] >= min_size:
            result[output == i + 1] = 255
    return result

def get_kernel_disk (imclose_para):
    if imclose_para == 5:
        return np.array(
[
   [0,   0,   1,   1,   1,   1,   1,   0,   0],
   [0,   1,   1,   1,   1,   1,   1,   1,   0],
   [1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1],
   [0,   1,   1,   1,   1,   1,   1,   1,   0],
   [0,   0,   1,   1,   1,   1,   1,   0,   0]
],dtype=np.uint8)

    elif imclose_para == 7:
        return np.array([
   [0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0],
   [0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1],
   [0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0],
   [0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0]
        ]
,dtype=np.uint8)

    elif imclose_para == 15:

        return np.array([
   [0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,  0,   0,   0,   0],
   [0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,  0,   0,   0,   0],
   [0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,  0,   0,   0,   0],
   [0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,  0,   0,   0,   0],
   [0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  0,   0,   0,   0],
   [0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   0,   0,   0],
   [0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   0,   0],
   [0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   0],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   1],
   [0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   1,   0],
   [0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   1,   0,   0],
   [0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  1,   0,   0,   0],
   [0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,  0,   0,   0,   0],
   [0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,  0,   0,   0,   0],
   [0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,  0,   0,   0,   0],
   [0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,  0,   0,   0,   0],
   [0,   0,   0,   0,   0,   0,   0,   0,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,  0,   0,   0,   0]
  ]
,dtype=np.uint8)
