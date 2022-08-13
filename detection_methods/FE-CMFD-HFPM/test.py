from process_image import Parameters
from feature_matching  import *
from post_processing import *
import matplotlib.pyplot as plt
from glob import glob
from tqdm.contrib.concurrent import process_map
from pathlib import Path
from PIL import Image
import cv2
import os

files = glob("/work/jcardenuto/translate-forgery-algs/testset/**/figure.png", recursive=True)

def run_method(fname):

    #img = cv2.imread(fname)
    #img_rgb = img[...,::-1]
    img_rgb = Image.open(fname)
    img_rgb = np.array(img_rgb)

    # set_parameters
    para = Parameters(img_rgb)
    para.plotimg = 0

    p1, p2 = feature_matching(fname, para)
    para.p1, para.p2 = p1.copy(),p2.copy()


    para.check_orienation = False
    ret = post_processing(img_rgb, para)

    map = ret[1]
    final_map = generated_object_match_map(map, p1, p2)


    vis_image = img_rgb.copy()
    if len(p1):
        for p_index in range(p1.shape[1]):
            point1 = (int(p1[0,p_index]), int(p1[1,p_index]))
            point2 = (int(p2[0,p_index]), int(p2[1,p_index]))

            vis_image = cv2.line(img=vis_image , pt1=point1 , pt2=point2, color=(0,255,0), thickness=2)
    else:
        vis_image = img_rgb


    fname = fname.split('testset/')[-1]
    fname = fname[:-4]
    os.makedirs(f'result/{fname}',exist_ok=True)

    np.save(f"result/{fname}_map.npy", map)
    np.save(f"result/{fname}_final_map.npy", final_map)
    plt.imsave(f'result/{fname}_final_map.png', final_map)
    plt.imsave(f'result/{fname}_matched_image.png', vis_image)


run_method("/work/jcardenuto/translate-forgery-algs/testset/compound/inter-panel/duplication/splicing/00996/figure.png")
#process_map(run_method, files[:100], chunksize = 10, max_workers = 20)

