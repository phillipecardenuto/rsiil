from re import L
from time import sleep
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
import cv2
from pywt import dwt2
from cyvlfeat import sift
from skimage.segmentation import slic, mark_boundaries
from skimage.measure import regionprops
from numpy.linalg import norm
from scipy.spatial.distance import pdist, squareform, cdist

def getCircleMask(B):
	circle=np.zeros((B,B));
	for i in range(B):
		for j in range(B):
			if norm([i-int((B+1)/2),j-int((B+1)/2)])<=B/2:
				circle[i,j]=1;
	return circle.astype("uint8")


def  map2RGB(map,RGBimage,color):
	red = RGBimage[:,:,0]
	green = RGBimage[:, :,1];
	blue = RGBimage[:,:,2];
	red[map] = color[0];
	green[map] = color[1];
	blue[map] = color[2]
	UltimateResult = RGBimage;
	UltimateResult[:,:,0] = red;
	UltimateResult[:,:,1] = green;
	UltimateResult[:,:,2] = blue
	return UltimateResult



def difference_angular(alpha1,alpha2):
	a1 = np.mod(alpha1,360);
	a2 = np.mod(alpha2,360);
	a1 = a1+(a1<0)*360;
	a2 = a2+(a2<0)*360;
	b1 = abs(a1-a2);
	b2 = min(a1,a2)+360-max(a1,a2);
	out = min(b1,b2);
	return out

def oversegmentation(filename):
	#@ Return Forgery Map , Match List , Locations

	color=[64,0,128];
	TR_P=2;
	TR_sim=15;

	RGBimage=np.array(Image.open(filename)).astype("uint8")
	M,N = RGBimage.shape[:2]
	if len(RGBimage.shape) > 2:
		# Discard alpha layer if exists:
		if RGBimage.shape[2] >3:
		     RGBimage = RGBimage[:3]
		grayimage = cv2.cvtColor(RGBimage.astype(np.uint8), cv2.COLOR_RGB2GRAY).astype(np.float64);
	else:
		grayimage = RGBimage

	Final_map1=np.zeros((M,N));
	Final_map2=np.zeros((M,N));
	Final_map3=np.zeros((M,N));

	CA1,(CB1,CC1,CD1)=dwt2(RGBimage,'haar');
	CA2,(CB2,CC2,CD2)=dwt2(CA1,'haar');
	CA3,(CB3,CC3,CD3)=dwt2(CA2,'haar');
	CA4,(CB4,CC4,CD4)=dwt2(CA3,'haar');
	E_LF=np.sum(abs(CA4));
	E_HF=np.sum(abs(CB1))+np.sum(abs(CB2))+np.sum(abs(CB3))+np.sum(abs(CB4))+ \
	    	np.sum(abs(CC1))+np.sum(abs(CC2))+np.sum(abs(CC3))+np.sum(abs(CC4))+ \
	    	np.sum(abs(CD1))+np.sum(abs(CD2))+np.sum(abs(CD3))+np.sum(abs(CD4));
	P_LF=E_LF/(E_LF+E_HF)*100;
	S=(P_LF>50)*np.sqrt(0.02*M*N)+(P_LF<=50)*np.sqrt(0.01*M*N);
	segments = slic(RGBimage,n_segments=S, start_label=1);
	num_segments= np.max(segments);


	"""
	if show==1:
		boundaries =  cv.Laplacian(segments, ddepth=-1, ksize=3)>0.1
		boundaries=repmat(boundaries,[1,1,3]);
		boundaries=imdilate(boundaries,ones(5));
		imp = RGBimage;
		imp(boundaries) = 0 ;
		figure;
		imshow(imp);
	end
	"""


	### SIFT
	Locations, Descriptors = sift.sift( grayimage.astype(np.float),
					     n_octaves = 5,
				            edge_thresh = 12,
					    peak_thresh = 0.1,
					    compute_descriptor=True)
	#Locations[:, [1,0]] = Locations[:, [0,1]]

	num_keypoint = Locations.shape[0];
	ratio=num_keypoint/(M*N);

	### Descriptor Normalization
	Descriptors = Descriptors.astype(np.float64);
	Descriptors = Descriptors / norm(Descriptors, axis=0, ord=2)

	### Matching
	### % Assign a segment number to each keypoint
	keypoint_segment= np.zeros(num_keypoint);
	for k in range(num_keypoint):
		keypoint_segment[k] = segments[ np.round(Locations[k,0]).astype(np.int),np.round(Locations[k,1]).astype(np.int)]

	### STEP1: Compare each 2 segments and create a matchList
	CC=[] #np.zeros(int(num_segments*(num_segments-1)/2));
	MatchList=[];
	num=-1;
	CC = np.zeros(int(num_segments*(num_segments-1)/2))
	for i in range(num_segments):
		# Superpixels labeled with the same id
		idx1=np.where(keypoint_segment==i)[0];
		if len(idx1)<2:
			num+=(num_segments-i-1);
			continue;
		D1=Descriptors[idx1,:];

		min1=np.sort(
				squareform(
					(pdist(D1))),axis=1);
		min1=min1[:,1] + np.finfo(float).eps # Sum eps to avoid division by zero

		for j in range(i+1,num_segments):
			num+=1;
			# Superpixels labeled with the same id
			idx2=np.where(keypoint_segment==j)[0];

			if len(idx2)<2:
				continue;
			D2=Descriptors[idx2, :];
			#min2=min(pdist(D2));
			#min3=min(min1,min2);
			feat_dist = cdist(D1,D2);  # %Distances between descriptors
			feat_dist, idx= np.sort(feat_dist, axis=1), np.argsort(feat_dist, axis=1); #%Sort rows
			idx=idx[:, 0]; # get all nearest kp (kpi <-> kpj)
			feat_dist=feat_dist[:,0]/min1; #feat_dist(:,2);%Ratio of 1st minima to 2nd
			k1=np.where(feat_dist<=1/TR_P)[0];

			new_pairs = 0
			if len(k1):
				# Insert point labeled as i and j that has the nearest features desc
				for n in k1:
					kp1, kp2 = idx1[n], idx2[idx[n]]#	new_pairs=np.hstack([idx1[k1] , idx2[idx[k1]]]);
					MatchList.append((kp1, kp2));
					new_pairs += 1
			CC[num] = new_pairs #   if len(k1) else 0

	num_match=len(MatchList);
	if num_match==0:
		return np.zeros_like(grayimage), [], []

	### %Show Match Keypoints
	"""
	if (show==1):
		figure;
		imshow(RGBimage);
		hold on;
		title('Match Keypoints');
		loc1_temp=Locations(MatchList(:,1),1:2);
		loc2_temp=Locations(MatchList(:,2),1:2);
		locs_temp=[loc1_temp;loc2_temp];
		plot(locs_temp(:,2),locs_temp(:,1),'mo');%Show points
		temp=[loc1_temp,loc2_temp,nan(size(loc2_temp,1),4)]';
		locs_temp=reshape(temp(:),4,[]);
		plot([locs_temp(2,:);locs_temp(4,:)],[locs_temp(1,:);locs_temp(3,:)],'b-');
		hold off;
	end
	"""
	# %STEP2: Calculate the TR_B
	CC_Sorted=np.sort(CC,axis=0);
	der_1st= np.gradient(CC_Sorted)# cv2.filter2D( CC_Sorted, kernel=np.array([[-1],[1]]), ddepth=-1,borderType=cv2.BORDER_REPLICATE);
	mean_der_1st=np.mean(der_1st);
	der_2nd= np.gradient(der_1st) #cv2.filter2D(CC_Sorted,kernel=[[-1],[2],[-1]], ddepth=-1, borderType=cv2.BORDER_REPLICATE);
	TR_B=min(CC_Sorted[der_2nd>mean_der_1st])+1 if len(CC_Sorted[der_2nd>mean_der_1st]) >0 else 1;#%%%%+2

	## %STEP3: Removing blocks with low matched points
	CC_Square=squareform(CC)>TR_B;

	# %STEP4: removing matched points such that CC_Square==0
	select=np.zeros(num_match);
	for k in range(num_match):
		s1=int(keypoint_segment[MatchList[k][0]]);
		s2=int(keypoint_segment[MatchList[k][1]]);
		select[k]=CC_Square[s1,s2];

	select = select !=0
	new_MatchList=np.array(MatchList)[select,:];
	num_match=new_MatchList.shape[0];
	if num_match==0:
		return np.zeros_like(grayimage), [], []

	#%% Forgery Region Extraction
	#%STEP1: new segmentation and assign segments to each matched keypoint
	superpixel_size=(max(M,N)>=1500)*10+10;
	if max(M,N) >=1500:
		n = 500
	else:
		n = 100
	superpixel = slic(RGBimage,n ,start_label=1);
	num_superpixel=np.max(superpixel);
	# % Assign a superpixel number to each matched keypoint
	Match_superpixel = np.zeros((num_match,2));#%SR
	superpixel_select = np.zeros(num_superpixel+1); # %Which superpixels have match and must processed
	for k in range(num_match):
		t1=new_MatchList[k,0];
		t2=new_MatchList[k,1];
		Match_superpixel[k,0] = superpixel[int(np.round(Locations[t1,0])),int(np.round(Locations[t1,1]))];
		Match_superpixel[k,1]=superpixel[int(np.round(Locations[t2,0])),int(np.round(Locations[t2,1]))];
		superpixel_select[int(Match_superpixel[k,0])]=1;
		superpixel_select[int(Match_superpixel[k,1])]=1;

	superpixel_select_list= np.where(superpixel_select)[0]; #%Convert logical vector to index
	superpixel_neighbor=np.zeros((num_superpixel+1,8));
	superpixel_center=np.zeros((num_superpixel+1,2));
	superpixel_color=np.zeros(num_superpixel+1);
	# %STEP2: find neighbors of matched segments
	# %computing the center of each superpixel
	regions = regionprops(superpixel, intensity_image=grayimage)
	for i, loc in enumerate(regions):
		y,x = loc.centroid
		superpixel_center[i,0] = y
		superpixel_center[i,1] = x
	"""
	for i in range(num_superpixel):
		map= superpixel==i;
		y,x=np.where(map);
		superpixel_center[i,0]=np.mean(y);
		superpixel_center[i,1]=np.mean(x);
	"""

	#%finding neighbors + create initial final_map

	for i in range(len(superpixel_select_list)):
		map= superpixel==superpixel_select_list[i];
		Final_map1[map]=1;

	    #%Finding neighbor superpixels
		center = superpixel_center[superpixel_select_list[i],:];
		kernel = cv2.getStructuringElement(cv2.MORPH_RECT,(7,7))
		dil_map = cv2.morphologyEx(map.astype("uint8"), cv2.MORPH_OPEN, kernel)
		y,x = np.where(dil_map-map>0);
		neighbors=np.zeros(num_superpixel+1);
		for j in range(len(x)):
			temp=superpixel[y[j],x[j]];
			neighbors[temp]=1;


		idx1 = np.where(neighbors)[0];
		neighbor_centers = superpixel_center[idx1,:];
		neighbor_angles = np.arctan2(neighbor_centers[:,0]-center[0],neighbor_centers[:,1]-center[1])*180/np.pi;
		for j in range(8):
			idx2=np.where(difference_angular(neighbor_angles,j*45-180)<=22)[0];
			if len(idx2)==0:
				continue;
			superpixel_neighbor[superpixel_select_list[i],j]=idx1[idx2[0]];

	# %Compute mean color of found neighbor superpixels
	temp= np.zeros(num_superpixel+1);
	temp[np.any(superpixel_neighbor,axis=1)] = 1;
	temp = np.where(temp)[0];
	for i in range(len(temp)):
		#TODO debug the next repeat line
		map=superpixel==temp[i]
		# RGB Image is actually RGB
		if len(RGBimage.shape) > 2:
		    map = np.tile(map, (3,1,1)).reshape(map.shape[0], map.shape[1],3)
		superpixel_color[temp[i]] = np.mean(RGBimage[map]);
	#%Compare neighbors using its colors
	Final_map2 = Final_map1.copy();
	for i in range(num_match):
		s1=int(Match_superpixel[i,0])
		s2=int(Match_superpixel[i,1])
		for j in range(8):
			n1=int(superpixel_neighbor[s1,j])
			n2=int(superpixel_neighbor[s2,j])
			if n1==0 or  n2==0:
				continue;

			if abs(superpixel_color[n1]-superpixel_color[s1])<=TR_sim and \
				 abs(superpixel_color[n2]-superpixel_color[s2])<=TR_sim:

				Final_map2[np.logical_or(superpixel==n1, superpixel==n2)] = 1
    #%STEP3: Morphological Operation

	kernel = getCircleMask(2*superpixel_size)
	Final_map2 = Final_map2.astype("uint8")
	Final_map3 = cv2.morphologyEx(Final_map2, cv2.MORPH_CLOSE, kernel)

	return Final_map3 , MatchList , Locations
