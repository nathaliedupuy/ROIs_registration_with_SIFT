# ROIs_registration_with_SIFT
Register the ROIs across multiple sessions for a given animal and a given field of view. 
The fields-of-view are aligned using Scale-Invariant Feature Transform (SIFT; David Lowe 1999).
ROIs matching is based on the method implemented in CaImAn (Giovannucci et al. eLife 2019).

## Requirements
- MATLAB and Python. The main code and GUI are in MATLAB, while some features are in Python (image filters, SIFT).
- See INSTALLATION_GUIDE.md for a step-by-step installation.

## Registration steps
	(0)	Prepare the dataset: template image and rois mask for each session.
		This part was developed to extract the data processed according to the rtools pipeline
		(average tiff images of all recordings and extract rois coordinates from imageJ zip file).
	(1)	Image preprocessing (enhance contrast and sharpness of image).
	(2)	Align images with SIFT (find local features and match them).
	(3)	Register the ROIs between sessions.
	(4)	GUI allows the user to edit the dataset with registration results.

## External licenced code 
- NoRMCorre: motion correction algorithm developed by the CaImAn team (Pnevmatikakis, 2016, https://github.com/flatironinstitute/NoRMCorre).
Used in step (0) to prep the template image of each session before registration.

- Hungarian : optimization algorithm (see https://en.wikipedia.org/wiki/Hungarian_algorithm). Used in step (3) to find an optimal matching between the ROIs.
The registration between two sessions is based on the registration implemented by the CaImAn team. 
For more info see their original function register_ROIs.m in their GitHub repo (https://github.com/flatironinstitute/CaImAn-MATLAB) and method in their paper (Giovannucci et al. eLife 2019 https://elifesciences.org/articles/38173).

- ReadImageJROI from Dylan Muir (https://github.com/DylanMuir/ReadImageJROI). Used in step (0) to extract Rois coordinates from ImageJ rois zip files.

- OpenCV image filters to enhance image.

- SIFT (library in OpenCV): scale-invariant feature transform is a powerful algorithm for local feature detection in computer vision developed by David Lowe in 1999. It was patented but still has free access for research.
