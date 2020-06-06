
%% README  
%
% This is the main script to register multiple sessions for a given animal
% and a given field of view. 
%
% The steps can be run independently:
% (0)   [NRLab] Prepare the dataset: template image and rois mask for each
%       session  (average recs tiff of and extract rois coordinates from
%       imageJ zip file)
% (1)   Prepro image (enhance contrast and sharpness of image)
% (2)   Align images with SIFT (find features and match them)
% (3)   Register the ROIs between sessions.
% (4)   The GUI allows the user to edit the dataset with registration results.
%
% _________________________________
%
%  External licenced code 
% _________________________________
%
% - NoRMCorre: motion correction algorithm developed by the CAIMAN team.
% (Pnevmatikakis, 2016)   https://github.com/flatironinstitute/NoRMCorre
% Used to prep the template image of each session before registration
%
% - Hungarian : optimization algorithm (see
% https://en.wikipedia.org/wiki/Hungarian_algorithm). Used in step (3) to
% find an optimal matching between the ROIs. The registration between two
% sessions is based on the registration implemented by the CaImAn team. For
% more info see their original function register_ROIs.m in their GitHub
% repo (https://github.com/flatironinstitute/CaImAn-MATLAB) and method in
% their paper (Giovannucci et al. eLife 2019
% https://elifesciences.org/articles/38173).

% - ReadImageJROI from Dylan Muir : Extract Rois coordinates ImageJ 
%       https://github.com/DylanMuir/ReadImageJROI
% 
% - SIFT (using implementation in OpenCV): scale-invariant feature
% transform is an algorithm for feature detection in computer vision
% developed by David Lowe in 1999. It was patented but still has free
% access for research.
%
%


%% PATH AND UNIQUE_ID FOR REGISTRATION

% Create UID for the registration project. This will be the name of the
% folder that contains the registration data for the given sessions.
registration_UID = ''; 

% List of sessions corresponding to the animal+fov to register
list_sessions = {'', ''};  % add as many as required

% Path to the parent folder for registration, which will contain the folder
% named 'registration_UID'
path_regdata='path\to\REGDATA';


%% (0) PREP DATASET (customised for NRLab)
%
% This code section generates:
%   - sesname_original.tif  
%           Tiff file of the original fov image - format32 float
%   - sesname_rois.mat
%           The rois coordinates and IDs - contains two variables named:
%           ROIs_XY and ROIsIDs

% Extract from working directory all specified sessions for a given
% animal/fov. For each session, motion correct (optional) then average
% across trials to create template image. Export ROIs coordinates from
% the imageJ zip files. The dataset is saved in the folder
% 'registration_UID'.

path_workingdir = ''; % Full path to working directory 
the_animal = '';
the_fov = []; % Leave empty if for some reason this name changes from session to session...
% ! Note that all recording tiffs that are present in a session folder will
% be loaded.

% Setup options
options.prepdata.export_rois = true; % Note 
options.prepdata.rois_2exclude=[1,2]; % Remove globals (background and glob npil)
options.prepdata.export_images = true; % Ignore all options below if false
options.prepdata.do_align_across = true; % set to false if already ran mc across in main pipeline BUT make sure the tiff corrected aross are in the main folder
% NOTE: it is recommended to re-correct across if you had previously done
% align across with sima.
% Optional, for motion correction if do_align_across=true :
options.prepdata.mc_max_shift = 100; % in pixels
options.prepdata.print_mc_msg = false; % only set to true if you wish to display the messages from motion correction

tic
prepare_data_nrlab(path_regdata, registration_UID, path_workingdir, the_animal, the_fov, list_sessions, options.prepdata);
toc


%% (1,2,3) REGISTRATION
%
% Each step can bu run independently.
% (1)   Prepro image (enhance contrast and sharpness of image)
% (2)   Align images with SIFT (find features and match them)
% (3)   Register the ROIs between sessions.

options.registration.root_path = '......\ROIs_registration_with_SIFT'; % Path to roi registration tool
options.registration.python_path = '....\python';  % Full path to python  (if you use virtual environment and started MATLAB from base)
options.registration.settingsfile = 'setting_python_regROIs.yaml';


% (1)   Prepro image (enhance contrast and sharpness of image)
% Requirements:
% --------------
%   - sesname_original.tif
%           Tiff file of the original fov image - format32 float and scaled
%           in range [0,255]

options.registration.prep_images = true; % If false, skip step (1)
options.registration.do_enhance_template = true; % Ignore if prep_images=false

% (2)   Align images with SIFT (find features and match them) and save
% results in subfolder   registration_UID/KEYPT_DATA
% Requirements:
% --------------
%   - sesname.tif
%           Tiff file of the fov image converted and enhanced (optional)
%           used for SIFT (sesname.tif) - format uint8 and scaled [0,255]
%           See step 1.

options.registration.align_fovs = true; % If false, skip step (2)

% >> See settingsfile to set options for image matching <<

% (3)   Register the ROIs between sessions.
% Requirements:
% --------------
%   - sesname.tif
%           Tiff file of the fov image converted and enhanced (optional)
%           used for SIFT (sesname.tif) - format uint8 and scaled [0,255]
%           See step 1.
%   - sesname_rois.mat
%           The rois coordinates and IDs - contains two variables named:
%           ROIs_XY and ROIsIDs
%   - SIFT keypoints in folder registration_UID/KEYPT_DATA, as generated
%           during step 2. See 'run_ROIsRegistration' for more info.


options.registration.register_rois = true; % If false, skip step (3)

% >> Options below can be ignored if register_rois=false <<
% Detect the ROIs that are potentially out of FOV:
options.registration.roireg.min_roi_area = 0.5; % Area_fov2 >= p * Area_fov1
% Detect bad keypoints using the inferred global transformation: recommended
options.registration.roireg.rmv_kp_outliers = true; 
options.registration.roireg.max_kpTest_dpix = 10; % distance inferred location vs actual
% Improve ROIs alignment with local transformation: recommended
options.registration.roireg.local_match = true; 
% Matching options (see CAIMAN)
options.registration.roireg.match_opt.dist_thr = 0.5; % see 'matchROIs_2ses.m'
% Change the ROIs contour into ellipses (improves overlap and thus matching)
options.registration.roireg.use_ellipses = true;
options.registration.roireg.ellipse_size_constrain=1.1; % the inferred ellipse cannot be more than x * original_axis_length



tic
[REG_RES, REG_DATA, KeyPoints] = run_ROIsRegistration(path_regdata, registration_UID, list_sessions, options.registration);
toc



 



%% GUI TO INSPECT AND EDIT RESULTS
%
% Modifies the results and save it into a new structure
% - Add/Delete pairs
% - Flag ROIs out of FOV
% - Delete ROIs
% 

updated_struct_name = 'REG_RES_updated';

GUI_regROIs(REG_DATA, REG_RES, updated_struct_name);


% Once you are done and happy with the results, save the last structure and
% you can delete the rest (including the registration folders). 

% If you wish to keep track of the registration process you could always
% save 'REG_DATA', 'REG_RES' and the final updated registration structure,
% and still deleted all other data (ie keypoints, tif images etc..) since
% all you need to run the GUI is in 'REG_DATA'.



%% VISUALISE KEYPOINTS 

iReg = 1;

figure(100)
clf
hold on
for k=1:size(KeyPoints{iReg}.kp{2},1)
    plot(KeyPoints{iReg}.kp{2}(k,1), KeyPoints{iReg}.kp{2}(k,2),'MarkerEdgeColor','none','Marker','o','MarkerFaceColor', [0.3,0.3,0.8], 'MarkerSize', 6)
    plot([KeyPoints{iReg}.kp{2}(k,1), KeyPoints{iReg}.kp{1}(k,1)],[KeyPoints{iReg}.kp{2}(k,2), KeyPoints{iReg}.kp{1}(k,2)], 'Color', [0.5,0.5,1], 'LineWidth',2)
end
for k=1:size(KeyPoints{iReg}.kp_discard{2},1)
        plot(KeyPoints{iReg}.kp_discard{2}(k,1), KeyPoints{iReg}.kp_discard{2}(k,2),'MarkerEdgeColor','none','Marker','o','MarkerFaceColor', [0.8,0.3,0.3], 'MarkerSize', 6)
        plot([KeyPoints{iReg}.kp_discard{2}(k,1), KeyPoints{iReg}.kp_discard{1}(k,1)],[KeyPoints{iReg}.kp_discard{2}(k,2), KeyPoints{iReg}.kp_discard{1}(k,2)], 'Color', [1,0.5,0.5], 'LineWidth',2)
end
axis ij
axis equal


