import os
import argparse
import sys
import yaml
import cv2
import numpy as np
from scipy.io import savemat
import custom_loadmatfiles as cloadmat


# NOTE: to ensure true division, make it float!! or use np.true_divide


###############################
# Filters
###############################

# Filter image - smooth and sharpen image
def filter_fov(img, kernel_sharp, d, colorspread, geospread):
	# Filters - smooth / sharp
	# Better to smooth/sharpen before contrast stretch!
	# Check: Unsharp masking
	# Used: Bilateral filter (lowpass) + sharp (highpass)

	'''NOTE for Kernel (wiki)
	Normalization -  division of each element in the kernel by the sum of all kernel elements,
	so that the sum of the elements of a normalized kernel is one.
	This will ensure the average pixel in the modified image is as bright as the average pixel in the original image
	'''
	img_filt = cv2.bilateralFilter(img,d,colorspread,geospread)
	img_filt_sharp = cv2.filter2D(img_filt, -1, kernel_sharp)
	return img_filt_sharp

'''
Bilateral filtering: paper Tomasi & Manduchi 1998
http://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Tomasi98.pdf

NICE link : https://people.csail.mit.edu/sparis/bf_course/slides08/03_definition_bf.pdf

https://docs.opencv.org/2.4/doc/tutorials/imgproc/gausian_median_blur_bilateral_filter/gausian_median_blur_bilateral_filter.html

Parameters:
d: The diameter of each pixel neighborhood.
sigma_Color: Standard deviation in the color space.
sigma_Space: Standard deviation in the coordinate space (in pixel terms)

- geometric spread sigma_Space (sigma_d in paper) - chosen based on the desired amount of
low-pass filtering (large blurs more, i.e. combines values from more
distant image locations)

- photometric spread sigma_Color (intensities, sigma_r in paper) is set to achieve the desired
amount of combination of pixel value (i.e. pixels with values much closer
to each other than sigma_Color are mixed together and values much more
distant than sigma_Color are not. IMPORTANT: sigma_Color depends on re-scale
intensities!! so double check if use normalized values in range [0,1].
'''

# Contrast stretch across image (sliding window)
def correctImgCS(img, binsize, numbins, apha_low=0, alpha_high=1):

    d1, d2,  = img.shape
    list_krow = np.linspace(0, d1-binsize['rows'], num=numbins['rows'], endpoint=True, dtype=int)
    list_kcol = np.linspace(0, d2-binsize['columns'], num=numbins['columns'], endpoint=True, dtype=int)
    img_patch_cs = np.zeros((d1, d2))
    counter_pix_cs = np.zeros((d1,d2))

    for k_r in range(0,numbins['rows']):
        r_stop = list_krow[k_r]+binsize['rows']
        for k_c in range(0,numbins['columns']):
            c_stop = list_kcol[k_c]+binsize['columns']
            img_patch = img[list_krow[k_r]:r_stop, list_kcol[k_c]:c_stop]
            img_patch_cs[list_krow[k_r]:r_stop,
                         list_kcol[k_c]:c_stop] = img_patch_cs[list_krow[k_r]:r_stop,list_kcol[k_c]:c_stop] + contraststretch255(img_patch,apha_low=apha_low,alpha_high=alpha_high)[0]
            counter_pix_cs[list_krow[k_r]:r_stop,
                           list_kcol[k_c]:c_stop] = counter_pix_cs[list_krow[k_r]:r_stop, list_kcol[k_c]:c_stop] + 1

    img_full_cs = np.true_divide(img_patch_cs,counter_pix_cs)

    return img_full_cs


###############################
# Various utitlities
###############################

def contraststretch255(img, **kwargs):
	if 'apha_low' in kwargs.keys() and kwargs['apha_low']>0:
		img_min = np.quantile(img,kwargs['apha_low'])
	else:
		img_min = np.amin(img)

	if 'alpha_high' in kwargs.keys() and kwargs['alpha_high']<1:
		img_max = np.quantile(img,kwargs['alpha_high'])
	else:
		img_max = np.amax(img)

	img = minmaxrescale255(img, img_min, img_max)
	return img, img_min, img_max

def minmaxrescale255(img, img_min, img_max):

    img = np.true_divide(img - img_min , img_max - img_min)
    img = np.minimum(np.maximum(img,0), 1)
    img = 255 *img
    return img

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

def repatchimg(img, crop_edges):

	crop_row = any(crop_edges['rows'])
	crop_col = any(crop_edges['columns'])

	if not crop_row and not crop_col: return img

	if crop_col:
		cat_left = np.zeros((img.shape[0],crop_edges['columns'][0]))
		cat_right = np.zeros((img.shape[0],crop_edges['columns'][1]))
		img = np.concatenate((cat_left,img,cat_right),axis=1)

	if crop_row:
		img = np.concatenate((np.zeros((crop_edges['rows'][0],img.shape[1])),img,np.zeros((crop_edges['rows'][1],img.shape[1]))),axis=0)
		#img_rows_out = [ img_full[0:crop_edges['rows'][0],:] , img_full[-crop_edges['rows'][1]:,:] ]
		#img_cols_out = [ img_full[crop_edges['rows'][0]:-crop_edges['rows'][1],0:crop_edges['columns'][0]] , img_full[crop_edges['rows'][0]:-crop_edges['rows'][1],-crop_edges['columns'][1]:] ]
		#img_cols_out = [ img_full[:,0:crop_edges['columns'][0]] , img_full[:,-crop_edges['columns'][1]:] ]

	return img


###############################
# Core functions
###############################

def run_correction(img, prepro_settings):

	if prepro_settings['do_crop']:
		crop_edges = prepro_settings['crop_edges']
		# Crop - keep the edges as 0-arrays that will re-pad in the end
		img = cropimg(img, crop_edges)

	# Global contrast stretch to bring to range [0,255] where all parameters were assessed
	params_cs_global = prepro_settings['params_cs_global']
	img, img_min, img_max = contraststretch255(img, apha_low=params_cs_global['apha_low'], alpha_high=params_cs_global['alpha_high'])

	if prepro_settings['do_smooth_sharpen']:
		params_smoosh = prepro_settings['params_smoosh']
		# Filter smooth/sharp
		img = filter_fov(img, params_smoosh['kernel_sharp'], params_smoosh['d'], params_smoosh['colorspread'], params_smoosh['geospread'])

	if prepro_settings['do_contrast_stretch']:
		params_cs = prepro_settings['params_cs']
		# Adaptative contrast stretch across image
		img = correctImgCS(img, params_cs['binsize'], params_cs['numbins'], apha_low=params_cs['apha_low'], alpha_high=params_cs['alpha_high'])

	# Final global rescaling - this one is necessary if applied filters
	if prepro_settings['do_smooth_sharpen'] or prepro_settings['do_contrast_stretch']:
		img = contraststretch255(img)[0]

	if prepro_settings['do_crop']:
		return repatchimg(img, crop_edges)
	else:
		return img
	



def enhance_mlt_imgs(prepro_settings, path_source, path_output, save_info):

	# Loop to find the tiff files
	list_tifffiles = [f for f in os.listdir(path_source) if os.path.isfile(os.path.join(path_source, f)) and f[-4:] == '.tif']

	print("\t Running for "+str(len(list_tifffiles))+" files")

	# Check folder exists
	if not os.path.exists(path_output):
		os.makedirs(path_output)

	for f_tiff in list_tifffiles:

		# Load image
		img_original = cv2.imread(os.path.join(path_source, f_tiff), -1)

		img = run_correction(np.float32(img_original), prepro_settings)

		if prepro_settings['convert_uint8']:
			# Round and convert -
			img = np.uint8(np.round(img))
			'''
			NOTE: this creates an image of smaller size, while preserving information (see tests in jupyter)
			Also, SIFT works with this data type
			'''
		else:
			img = np.float32(img)

		# Save image
		cv2.imwrite(os.path.join(path_output, f_tiff), img)

	if save_info:
		# Save the parameters used to modify image
		#the_folder_an_output = os.path.join(path,test_output_name+'.mat')
		savemat(os.path.join(path_output,'PREPRO_INFO.mat'), prepro_settings)


def enhance_one_img(prepro_settings, path_source, name_source, path_output, name_output, save_info):

	# Load image
	img_original = cv2.imread(os.path.join(path_source,name_source+'.tif'), -1)
	img = run_correction(np.float32(img_original), prepro_settings)

	if prepro_settings['convert_uint8']:
		# Round and convert -
		img = np.uint8(np.round(img))
		'''
		NOTE: this creates an image of smaller size, while preserving information (see tests in jupyter)
		Also, SIFT works with this data type
		'''
	else:
		img = np.float32(img)
	
	
	

	# Save image
	if not path_output:
		path_output = path_source
		# It has to have a name output since it is in the same folder!
	else:
		if not name_output:
			name_output = name_source

	cv2.imwrite(os.path.join(path_output,name_output+'.tif'), img)

	if save_info:
		savemat(os.path.join(path_output, name_output+'_PREPRO_INFO.mat'), prepro_settings)



def main(settingsfile, path_source, name_source, path_output, name_output, save_info):

	if not settingsfile:
		# default parameters
		prepro_settings = {
			'convert_uint8':False,
			'do_crop':False,
			'do_smooth_sharpen':False,
			'do_contrast_stretch':True,
			'params_cs':{'apha_low':0, 'alpha_high':0.99, 'binsize':{'rows': 50, 'columns': 50}, 'numbins':{'rows': 200, 'columns': 200}},
			'params_cs_global':{'apha_low':0, 'alpha_high':1}}
		
	else:
		# user-defined settings
		# Open the parameters file
		with open(settingsfile, 'r') as hf:
			gbl_settings = yaml.safe_load(hf)

		# We just need the pre-pro settings
		prepro_settings = {
			'convert_uint8':gbl_settings['image_prepro']['convert_uint8'],
			'do_crop':gbl_settings['image_prepro']['do_crop'],
			'do_smooth_sharpen':gbl_settings['image_prepro']['do_smooth_sharpen'],
			'do_contrast_stretch':gbl_settings['image_prepro']['do_contrast_stretch'],
			'params_cs_global':gbl_settings['image_prepro']['params_cs_global']}

		# Add optional params in dictionary

		if prepro_settings['do_crop']:
			prepro_settings['crop_edges'] = gbl_settings['image_prepro']['crop_edges']

		if prepro_settings['do_contrast_stretch']:
			prepro_settings['params_cs'] = gbl_settings['image_prepro']['params_cs']
		
	# Check what function to call
	if not name_source:
		# Enhance FOVs images - read files in folder (recording of trials)
		#print("Enhance all images in folder "+path_source)

		if prepro_settings['do_smooth_sharpen']:
			# Setup the parameters for filters - add to gbl_settings
			# Basic sharpen kernel
			kernel_sharp = np.array([[-1,-1,-1], [-1,9,-1], [-1,-1,-1]], dtype=np.float32)
			prepro_settings['params_smoosh'] = {'kernel_sharp':kernel_sharp, 'd':5, 'colorspread': 75, 'geospread': 3}

		enhance_mlt_imgs(prepro_settings, path_source, path_output, save_info)

	else:
		#print("Enhance one image "+name_source)

		if prepro_settings['do_smooth_sharpen']:
			# Setup the parameters for filters - add to gbl_settings
			# Basic sharpen kernel
			kernel_sharp = np.array([[-1,-1,-1], [-1,9,-1], [-1,-1,-1]], dtype=np.float32)
			# May want to use Sharper for template:
			#kernel_sharp = np.array([[-1,-4,-1], [-4,21,-4], [-1,-4,-1]], dtype=np.float32)
			#kernel_sharp = np.array([[-1,-1,-1], [-1,9,-1], [-1,-1,-1]], dtype=np.float32)
			prepro_settings['params_smoosh'] = {'kernel_sharp':kernel_sharp, 'd':5, 'colorspread': 75, 'geospread': 3}

		enhance_one_img(prepro_settings, path_source, name_source, path_output, name_output, save_info)



if __name__ == '__main__':


    # Parse command line inputs
	#print 'Argument List:', str(sys.argv)
	parser = argparse.ArgumentParser(description='Enhance an fov image')
	parser.add_argument('settingsfile', nargs='?', default="")
	parser.add_argument('path_source', nargs='?', default="")
	parser.add_argument('name_source', nargs='?', default="")
	parser.add_argument('path_output', nargs='?', default="")
	parser.add_argument('name_output', nargs='?', default="")
	parser.add_argument('save_info', nargs='?', default=False)
	args = parser.parse_args()

	if args.save_info:
		if args.save_info=='1' or args.save_info=='True' or args.save_info=='true':
			args.save_info=True
		else:
			args.save_info=False

	# Run main
	main(args.settingsfile, args.path_source, args.name_source, args.path_output, args.name_output, args.save_info)



# NOT USED, BUT USEFUL!
'''
def check_whitespaces_in_path(*paths):
	no_wsp_paths = []
	for thepath in paths:
		if ' ' in thepath:
			no_wsp_paths.append(thepath.replace(' ','\ '))
		else:
			no_wsp_paths.append(thepath)
	return no_wsp_paths
'''
