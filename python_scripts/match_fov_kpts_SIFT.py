import os
import argparse
import sys
import yaml
import cv2
import numpy as np
from scipy.stats import pearsonr
from scipy.io import savemat
import custom_loadmatfiles as cloadmat


# Attention! images have to be [0,255] uint8 to work with SIFT
# Also assume images are of similar sizes (i.e. coordinate bases are similar, few pix differences)



###############################
# Prep image
###############################


def cropimg(img, crop_edges):

	crop_row = any(crop_edges['rows'])
	crop_col = any(crop_edges['columns'])

	if not crop_row and not crop_col: return img

	if crop_row and crop_col:
		img = img[crop_edges['rows'][0]:-crop_edges['rows'][1],crop_edges['columns'][0]:-crop_edges['columns'][1]]
	elif not crop_row and crop_col:
		img = img[:,crop_edges['columns'][0]:-crop_edges['columns'][1]]
	else:
		img = img[crop_edges['rows'][0]:-crop_edges['rows'][1],:]

	return img


def prepimg(prepro_settings, path_regdata, experiment_ID):

	# Load the settings to eventually crop image

	# Load prepo info for template if it exists
	if os.path.exists(os.path.join(path_regdata,experiment_ID+'_PREPRO_INFO.mat')):
		info_tpl = cloadmat.load_matlab_data(os.path.join(path_regdata,experiment_ID+'_PREPRO_INFO.mat'))

		# Check if template was cropped
		# Check if need to overwrite cropping settings
		if info_tpl['do_crop']:
			if prepro_settings['do_crop']:
				warn_crop = False
				for key in prepro_settings['crop_edges']:
					for k in [0,1]:
						precrop=info_tpl['crop_edges'][key][k]
						if not warn_crop and precrop>prepro_settings['crop_edges'][key][k]:
							print("WARNING - TEMPLATE CROPPING WAS LARGER - OVERWRITE CROP", flush=True)
							warn_crop=True
						prepro_settings['crop_edges'][key][k] = max(precrop, prepro_settings['crop_edges'][key][k])

			else:
				print("WARNING - TEMPLATE WAS CROPPED - FORCE CROP", flush=True)
				prepro_settings['do_crop'] = True
				prepro_settings['crop_edges'] = info_tpl['crop_edges']


	# Load image
	img = cv2.imread(os.path.join(path_regdata,experiment_ID+'.tif'), -1)

	# Get size
	originalsize=img.shape

	# Check cropping
	if prepro_settings['do_crop']:
		# Crop - keep the edges as 0-arrays that will re-pad in the end
		img = cropimg(img, prepro_settings['crop_edges'])
		prepro_settings['shift_xy'] = (prepro_settings['crop_edges']['columns'][0], prepro_settings['crop_edges']['rows'][0]) # crop left # crop top
	else:
		prepro_settings['shift_xy'] = (0,0) # crop left # crop top
		# REMEMBER X/Y -> IMAGE COL/ROW AND Y = flip up-down ROW

	return img, originalsize, prepro_settings


###############################
# Various utitlities
###############################

def get_match_KPxy(m, kp1, kp2):

	kp1_xy=(np.float32(kp1[m.queryIdx].pt[0]), np.float32(kp1[m.queryIdx].pt[1]))
	kp2_xy=(np.float32(kp2[m.trainIdx].pt[0]), np.float32(kp2[m.trainIdx].pt[1]))

	# Computed orientation of the keypoint (-1 if not applicable);
	# in [0,360) degrees and measured relative to image coordinate system, ie in clockwise.
	kp1_xy=kp1_xy+(np.float32(kp1[m.queryIdx].angle),) # add a singleton to be able to concat with tuple
	kp2_xy=kp2_xy+(np.float32(kp2[m.trainIdx].angle),)

	return kp1_xy, kp2_xy


def compute_match_KPdistance(kp1_xy,kp2_xy,shift_xy_1,shift_xy_2):

	dpix_xy = np.sqrt( np.float32(kp1_xy[0]+shift_xy_1[0] - kp2_xy[0]-shift_xy_2[0])**2 +
					   np.float32(kp1_xy[1]+shift_xy_1[1] - kp2_xy[1]-shift_xy_2[1])**2 )

	return dpix_xy

def compute_match_KPdxdy(kp1_xy,kp2_xy,shift_xy_1,shift_xy_2):

	dx = np.float32(kp1_xy[0]+shift_xy_1[0] - kp2_xy[0]-shift_xy_2[0])
	dy = np.float32(kp1_xy[1]+shift_xy_1[1] - kp2_xy[1]-shift_xy_2[1])

	return dx, dy


def MahalanobisDist(data, mu, covmat_inv):
	md = []
	for i in range(data.shape[1]):
		z = data[:,i] - mu
		md.append( np.float32( np.sqrt(np.dot(np.dot(z,covmat_inv),z)) ) )

	return md


def reset_KP_origcoord(kp_xy,shift_xy):

	if any(shift_xy):
		return ( np.float32(kp_xy[0]+shift_xy[0]) , np.float32(kp_xy[1]+shift_xy[1]), np.float32(kp_xy[2]) )
	else:
		return kp_xy


##################################
# Functions for matches selection
##################################

def getSquarePatch(x,y,siz,winpix,pix_shift_edg):

	# Transform coordinates into pixels
	x, y = [np.int(np.round(u)) for u in [x,y]]

	# Now get the rows and columns of patch
	rows_pix=[y-winpix, y+winpix+1]
	cols_pix=[x-winpix, x+winpix+1]
	create_patch = rows_pix[0]>=-pix_shift_edg and rows_pix[1]<=siz[0]+pix_shift_edg and cols_pix[0]>=-pix_shift_edg and cols_pix[1]<=siz[1]+pix_shift_edg
	if create_patch:

		if rows_pix[0]<0:
			shift_c_row=-rows_pix[0]
		elif rows_pix[1]>=siz[0]:
			shift_c_row=siz[0]-rows_pix[1]
		else:
			shift_c_row=0

		if cols_pix[0]<0:
			shift_c_col=-cols_pix[0]
		elif cols_pix[1]>=siz[1]:
			shift_c_col=siz[1]-cols_pix[1]
		else:
			shift_c_col=0

		patch = {'rows':rows_pix, 'columns':cols_pix}
		patch_shift = {'rows':shift_c_row, 'columns':shift_c_col}

	else:
		patch = {}
		patch_shift = {}

	return create_patch, patch, patch_shift


def getShiftPatch(patch1, patch2, patch1_shift, patch2_shift):

	temp_r=[patch1_shift['rows'],patch2_shift['rows']]
	temp_absr=[abs(patch1_shift['rows']),abs(patch2_shift['rows'])]
	temp_absr_max=max(temp_absr)

	shiftpatch_row0=temp_r[temp_absr.index(temp_absr_max)]

	patch1['rows']=[r+shiftpatch_row0 for r in patch1['rows']]
	patch2['rows']=[r+shiftpatch_row0 for r in patch2['rows']]

	temp_c=[patch1_shift['columns'],patch2_shift['columns']]
	temp_absc=[abs(patch1_shift['columns']),abs(patch2_shift['columns'])]
	temp_absc_max=max(temp_absc)

	shiftpatch_col0=temp_c[temp_absc.index(temp_absc_max)]

	patch1['columns']=[c+shiftpatch_col0 for c in patch1['columns']]
	patch2['columns']=[c+shiftpatch_col0 for c in patch2['columns']]

	return patch1, patch2, shiftpatch_row0, shiftpatch_col0




def getPatchSimilarity(img1, img2, kp1_xy, kp2_xy,
					   winpix=20, pix_shift_edg=0, min_corr=0.3):

	# *Inputs*
	# - winpix : square patch around kpt coordinates - patch size1=2*winpix+1 (always odd)
	# - pix_shift_edg : in case kpt is close to the edge, we shift the patch (! will no longer be centered on the kpts)

	check_edg=pix_shift_edg>0

	siz1=img1.shape
	siz2=img2.shape

	# Check if can create patches arround the keypoints - Note winpix is added on both sides of the kpt
	create_patch1, patch1, patch1_shift = getSquarePatch(kp1_xy[0],kp1_xy[1],siz1,winpix,pix_shift_edg)
	create_patch2, patch2, patch2_shift = getSquarePatch(kp2_xy[0],kp2_xy[1],siz2,winpix,pix_shift_edg)

	if not create_patch1 or not create_patch2 :
		return False, {}

	if check_edg:
		# Check if has an edge patch
		is_edg_patch = patch1_shift['rows']!=0 or patch1_shift['columns']!=0 or patch2_shift['rows']!=0 or patch2_shift['columns']!=0
		if is_edg_patch:
			# Shift the patches in case had to shift one
			patch1, patch2, shiftpatch_row0, shiftpatch_col0 = getShiftPatch(patch1, patch2, patch1_shift, patch2_shift)
		else:
			shiftpatch_row0, shiftpatch_col0 = [0, 0]

	else:
		is_edg_patch, shiftpatch_row0, shiftpatch_col0 = [False, 0, 0]

	# Square patches -> flatten
	sq_patch1=np.ndarray.flatten(np.float32(img1[patch1['rows'][0]:patch1['rows'][1],patch1['columns'][0]:patch1['columns'][1]]))
	sq_patch2=np.ndarray.flatten(np.float32(img2[patch2['rows'][0]:patch2['rows'][1],patch2['columns'][0]:patch2['columns'][1]]))

	# NOTE: with reshape, dot product does not work

	# Check initial correlation of the patch

	patch_sim=np.float32(pearsonr(sq_patch1,sq_patch2)[0])

	# Check if correlation is above minimum threshold
	if patch_sim<min_corr :
		# If not, Discard
		return False, {}
	else:
		# If it is already high, Automatically keep without adjusting
		return True, {'is_edg_patch':is_edg_patch, 'pcorr':patch_sim}

def remove_matches_outliers_MahalD(matches, high_corr, shift_xy_1, shift_xy_2, thresh_outlier=3):

	# Get the displacements of kpts between session 1 and 2 and compute Mahalanobis distance
	# >> Use the matches with good correlation to get the expected mean and covariance of MD

	# Find outliers using Mahalanobis distance
	'''see: https://stat.ethz.ch/education/semesters/ss2012/ams/slides/v2.2.pdf'''

	list_dx, list_dy=[],[]
	is_hc=[]

	# Loop over matches
	for m in matches:

		# Get keypoints coordinates and displacement
		if 'kp2_xy_adj' in m:
			dx, dy = compute_match_KPdxdy(m['kp1_xy'],m['kp2_xy_adj'],shift_xy_1,shift_xy_2)
			# Check correlation
			is_hc.append(m['pcorr_kp2adj']>high_corr)
		else:
			dx, dy = compute_match_KPdxdy(m['kp1_xy'],m['kp2_xy'],shift_xy_1,shift_xy_2)
			# Check correlation
			is_hc.append(m['pcorr']>high_corr)


		list_dx.append(dx)
		list_dy.append(dy)

	arr_dxdy = np.array([list_dx, list_dy]) # this will be of size [2 x numpts]
	arr_dxdy_high = np.array([[dx for dx,b in zip(list_dx,is_hc) if b], [dy for dy,b in zip(list_dy,is_hc) if b]]) # this will be of size [2 x numpts]

	dxdy_m=np.mean(arr_dxdy_high, axis=1)
	dxdy_cov=np.cov(arr_dxdy_high)
	''' Note for numpy cov(): 2-D array containing multiple variables and observations.
	Each row of m represents a variable, and each column a single observation of all those variables.'''
	# NOTE: could use the normal inv because we have a square mat, but just in case...
	dxdy_cov_inv=np.linalg.pinv(dxdy_cov)

	MD=MahalanobisDist(arr_dxdy, dxdy_m, dxdy_cov_inv)

	# Add it to the main dict
	[m.update(MahalD=m_md) for m,m_md in zip(matches, MD)]

	# Get the expected MD and std, using again the matches with high correlation
	list_MD_hc=[m_md for m_md, b in zip(MD,is_hc) if b]
	MD_hc_std=np.std(list_MD_hc)
	MD_hc_m=np.mean(list_MD_hc)

	return [m for m in matches if np.absolute(m['MahalD']-MD_hc_m) < thresh_outlier*MD_hc_std] # division ok because all are float



###############################
# Core functions
###############################



def run_match_2fovs(KEYPT_DATA, prepro_settings, path_regdata, expID_2align, sift, base_img, base_kp, base_des):

	# Prep the image to align
	alg_img, alg_originalsize, alg_prepro_settings = prepimg(prepro_settings, path_regdata, expID_2align)

	# Find the keypoints and descriptors with SIFT for the new image to align
	alg_kp, alg_des = sift.detectAndCompute(alg_img,None)

	# Update the dictionary
	KEYPT_DATA['info_alignedImg'] = {'name':expID_2align, 'size':alg_originalsize, 'prepro':alg_prepro_settings, 'SIFT_kpts': len(alg_kp)}

	base_shift_xy = KEYPT_DATA['info_baseImg']['prepro']['shift_xy']
	alg_shift_xy = alg_prepro_settings['shift_xy']

	# Get parameters to eval matches
	params_match = KEYPT_DATA['params_match']

	# Match the keypoints / features
	# Initiate BruteForce matcher with default params
	bf = cv2.BFMatcher()

	# Match descriptors.
	matches = bf.match(base_des, alg_des)

	# Only keep matches where kpt2 is within a set distance from kpt1
	# NOTE that it was recommended to select kp matches by comparing the distance
	# (i.e. euclidian distance between descriptors)
	# of the closest neighbor to that of the second-closest  neighbor (aka ratio test)
	# But we lose too many kp because they are likely similar with this type of images
	# ---> So instead we can apply other thresholds by using the fact that we expect a moderate deformation
	good_matches=[]
	for m in matches:

		kp1_xy, kp2_xy = get_match_KPxy(m, base_kp, alg_kp)

		dpix_xy=compute_match_KPdistance(kp1_xy, kp2_xy, base_shift_xy, alg_shift_xy)

		if dpix_xy < params_match['max_dxy_kpts'] :

			# Now  Inspect match : compare mini patches similarity (simple correlation)
			# Create a small patch around the kpts and compute similarity.
			# If similarity above min threshold but below good threshold, re-assess the kpt-match alignment
			# -> Transform patch into a round patch, and find the transformation (translating/rotating) of kpt2- patch that leads to best similarity.
			# Output the new coordinates+rotation, and similarity scores.

			keep_match, res = getPatchSimilarity(base_img, alg_img, kp1_xy, kp2_xy,
				params_match['patch_halfsiz'], params_match['patch_edg_shift'], params_match['patch_min_corr'])

			if keep_match:
				# Save the match!
				res['kp1_xy']=kp1_xy
				res['kp2_xy']=kp2_xy

				good_matches.append(res)

	# Save number of good matches found
	KEYPT_DATA['num_goodcorr_matches']=len(good_matches)

	# Final check, Mahalanobis distance to remove eventual outliers
	KEYPT_DATA['matches']=remove_matches_outliers_MahalD(good_matches, params_match['MD_high_corr'], base_shift_xy, alg_shift_xy, params_match['MD_outlier'])


	#print(KEYPT_DATA['matches'][0])
	# Reset KeyPoints in original image coordinate system
	[m.update(kp1_xy=reset_KP_origcoord(m['kp1_xy'],base_shift_xy)) for m in KEYPT_DATA['matches']]
	[m.update(kp2_xy=reset_KP_origcoord(m['kp2_xy'],alg_shift_xy)) for m in KEYPT_DATA['matches']]
	[m.update(kp2_xy_adj=reset_KP_origcoord(m['kp2_xy_adj'],alg_shift_xy)) for m in KEYPT_DATA['matches'] if 'kp2_xy_adj' in m]

	#print(KEYPT_DATA['matches'][0])

	# Save / export data
	expID_base = KEYPT_DATA['info_baseImg']['name']
	savemat(os.path.join(path_regdata,'KEYPT_DATA','KEYPT_DATA_'+expID_2align+'_to_'+expID_base+'.mat'), KEYPT_DATA)




def main(settingsfile, path_regdata, expID_base, expIDs_2align):

	# Open the parameters file
	with open(settingsfile, 'r') as hf:
		gbl_settings = yaml.safe_load(hf)


	# Initiate SIFT detector
	sift = cv2.xfeatures2d.SIFT_create()

	# We need the keypoints matching settings
	kpts_settings = gbl_settings['kpts_matching']

	print("Template session : " + expID_base, flush=True)

	# Prep the base image and info --
	base_img, base_originalsize, base_prepro_settings = prepimg(kpts_settings['prepro'], path_regdata, expID_base)

	# Find the keypoints and descriptors with SIFT for the base image (will be used to align the other images)
	base_kp, base_des = sift.detectAndCompute(base_img,None)

	# Prep the default output dict
	KEYPT_DATA = {'params_match':kpts_settings['params_match'], 'info_baseImg': {'name':expID_base, 'size':base_originalsize, 'prepro':base_prepro_settings, 'SIFT_kpts': len(base_kp)}}

	# Check folder exists
	if not os.path.exists(os.path.join(path_regdata,'KEYPT_DATA')):
		os.makedirs(os.path.join(path_regdata,'KEYPT_DATA'))

	for expid in expIDs_2align:

		print("Align session " + expid + " to " +expID_base, flush=True)

		run_match_2fovs(KEYPT_DATA, kpts_settings['prepro'], path_regdata, expid, sift, base_img, base_kp, base_des)



if __name__ == '__main__':

	# Parse command line inputs
	#print 'Argument List:', str(sys.argv)
	parser = argparse.ArgumentParser(description='Extract and match features with SIFT')
	parser.add_argument('settingsfile', nargs='?', default="")
	parser.add_argument('path_regdata', nargs='?', default="")
	parser.add_argument('expID_base', nargs='?', default="")
	parser.add_argument('expIDs_2align', nargs='+')
	args = parser.parse_args()

	# Run main
	main(args.settingsfile, args.path_regdata, args.expID_base, args.expIDs_2align)
