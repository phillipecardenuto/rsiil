from skimage import measure
import numpy as np
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