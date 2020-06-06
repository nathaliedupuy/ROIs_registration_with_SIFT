function [REG_RES, REG_DATA_1_match, REG_DATA_2_match, KeyPoints] = registerROIs_2sessions(path_regdata, session_1, REG_DATA_1, session_2, REG_DATA_2, options)
% Register ROIs between 2 sessions

% Potential TO_DO:
% Do a local alignement using matched ROIs and KP
% homography transformation

% REG_DATA structure with fields
%         'imgtiff'
%         'imgsiz' 
%         'ROIs_XY' % Rois coordinates are [Npts x 2(x,y)]
%         'ROIsIDs'
    

% -------------------------------------------------------------------------
% GET TRANSFORMATON FROM KEYPOINTS
% -------------------------------------------------------------------------

% Load KP data
KEYPT_DATA = load(fullfile(path_regdata, 'KEYPT_DATA', ['KEYPT_DATA_',session_2,'_to_',session_1,'.mat']));

% Get the keypoints coordinates
KeyPoints.kp = cell(1,2);
KeyPoints.kp{1}=cellfun(@(m)double(m.kp1_xy(1:2)), KEYPT_DATA.matches, 'UniformOutput', false); KeyPoints.kp{1}=cat(1,KeyPoints.kp{1}{:}); % [nPts x 2]
KeyPoints.kp{2}=cellfun(@(m)double(m.kp2_xy(1:2)), KEYPT_DATA.matches, 'UniformOutput', false); KeyPoints.kp{2}=cat(1,KeyPoints.kp{2}{:}); % [nPts x 2]

% Get Global affine transformation KP1 = M * KP2 + T
[M, T] = locfun_get_transformation_KP(KeyPoints.kp{1},KeyPoints.kp{2});
% size T: [2 x 1]

kp1_discard = []; kp2_discard = [];

% Check if we have outliers kpts matches by looking at the sse of their
% reconstruction using the transformation estimated above

if options.rmv_kp_outliers
    test_kp1 = KeyPoints.kp{2} * M' + T(:,ones(size(KeyPoints.kp{1},1),1))';
    se = sqrt(sum((KeyPoints.kp{1}-test_kp1).^2,2));
    good_kp = se < options.max_kpTest_dpix;
    
    % If we discarded kpts, re-calculate the global transformation
    if sum(good_kp)<size(KeyPoints.kp{1},1)
        warning(['Discarding ', num2str(sum(~good_kp)), ' KP match out of ', num2str(size(KeyPoints.kp{1},1))])
        kp1_discard = KeyPoints.kp{1}(~good_kp,:);
        kp2_discard = KeyPoints.kp{2}(~good_kp,:);
        KeyPoints.kp{1} = KeyPoints.kp{1}(good_kp,:);
        KeyPoints.kp{2} = KeyPoints.kp{2}(good_kp,:);
        [M, T] = locfun_get_transformation_KP(KeyPoints.kp{1},KeyPoints.kp{2});

    end
end
KeyPoints.kp_discard = {kp1_discard, kp2_discard};


% -------------------------------------------------------------------------
% TRANSFORM DATA - ROIS COORDINATES AND IMAGES
% -------------------------------------------------------------------------

% Transfo ROIs coordinates from session 2 to coord system of session 1
[REG_DATA_2_match.ROIs_XY.aligned] = locfun_transform_roi_coord(REG_DATA_2.ROIs_XY,M,T,'2_to_1');
% Rois coordinates are [Npts x 2(x,y)]

% Also do the inverse transformation inv(M)*(KP1 - T) = M \ (KP1 - T) = KP2
% so we can check the ROIs that are out of FOV (e.g. >50% area out)
[REG_DATA_1_match.ROIs_XY.aligned] = locfun_transform_roi_coord(REG_DATA_1.ROIs_XY,M,T,'1_to_2');



% Save a copy of rois coordinates. These data will be modifed later. See
% last section for explanation (SYNC part).
% - Original rois coordinates in the shared coordinate system (SYNC)
REG_DATA_1_match.ROIs_XY.original_sync = REG_DATA_1.ROIs_XY;
REG_DATA_2_match.ROIs_XY.original_sync = REG_DATA_2.ROIs_XY;
% - Aligned rois coordinates using the global transform only (SYNC)
REG_DATA_1_match.ROIs_XY.globalT_sync = REG_DATA_1_match.ROIs_XY.aligned;
REG_DATA_2_match.ROIs_XY.globalT_sync = REG_DATA_2_match.ROIs_XY.aligned;


% Get valid FOV: get the coordinates mapped to the other FOV and transform
% to binary masks; flag the ROIs whose area shrinks too much... i.e. are
% likely out of FOV. Remember we need the ROIs that are likely common to
% both fields of view!

[is_ok_bothfovs_1] = locfun_flag_outOfFov_ROIs(REG_DATA_1.ROIs_XY, REG_DATA_1_match.ROIs_XY.aligned, REG_DATA_2.imgsize, KEYPT_DATA.info_alignedImg.prepro, options.min_roi_area);
[is_ok_bothfovs_2] = locfun_flag_outOfFov_ROIs(REG_DATA_2.ROIs_XY, REG_DATA_2_match.ROIs_XY.aligned, REG_DATA_1.imgsize, KEYPT_DATA.info_baseImg.prepro, options.min_roi_area);

% Get the corners of the original and aligned images
[REG_DATA_1_match.corners_XY.original_sync, REG_DATA_1_match.corners_XY.aligned] = ...
    locfun_transform_imgcorners(M, T, REG_DATA_1.corners_XY, '1_to_2'); % [top left, bottom left, bottom right, top right]
[REG_DATA_2_match.corners_XY.original_sync, REG_DATA_2_match.corners_XY.aligned] = ...
    locfun_transform_imgcorners(M, T, REG_DATA_2.corners_XY, '2_to_1'); % [top left, bottom left, bottom right, top right]


% Tranform image s2 to align to s1
A = [[M',[0;0]] ; [T',1]];
tform = affine2d(A);
REG_DATA_2_match.imgtiff.aligned_sync = imwarp(REG_DATA_2.imgtiff,tform);

% Tranform image s1 to align to s2
A = [[inv(M)',[0;0]] ; [-(M\T)',1]];
tform = affine2d(A);
REG_DATA_1_match.imgtiff.aligned_sync = imwarp(REG_DATA_1.imgtiff,tform);

% Save a copy of images. These data will be modifed later. See
% last section for explanation (SYNC part).
REG_DATA_1_match.imgtiff.original_sync = REG_DATA_1.imgtiff;
REG_DATA_2_match.imgtiff.original_sync = REG_DATA_2.imgtiff;


% -------------------------------------------------------------------------
% MATCH ROIs IDs - CHECK IF REFINE ALIGNMENT WITH LOCAL TRANSORMATIONS
% -------------------------------------------------------------------------


if options.local_match
    % Refine alignment
    do_local_align();
end
% Match ROIs IDs - Note that we use the image 1 as template to project all
% ROIs, excluding the ones that are out of one of the FOVs.
[REG_RES] = matchROIs_2ses({REG_DATA_1.ROIs_XY(is_ok_bothfovs_1), REG_DATA_2_match.ROIs_XY.aligned(is_ok_bothfovs_2)}, ...
    {find(is_ok_bothfovs_1), find(is_ok_bothfovs_2)}, REG_DATA_1.imgsize, options.match_opt);

if options.use_ellipses
    [REG_RES_ell] = do_match_ellipse();
    
    % Merge the results
    [REG_RES] = locfun_merge_REGRES(REG_RES, REG_RES_ell);

end


    function do_local_align()
        
        
        radius_npix = 150;
        min_num_refpt = 8; % At least 3 matches are needed to provide a solution
        
        
        centroids_1 = cellfun(@(xy) mean(xy,1), REG_DATA_1.ROIs_XY, 'UniformOutput', false);
        centroids_2 = cellfun(@(xy) mean(xy,1), REG_DATA_2.ROIs_XY, 'UniformOutput', false);
        
        % TO DO :  ROIS =>
        % Get centroids for aligned data - this is the one that should be improved
        % using correlation of minipatch
        % centroids_2_aligned = cellfun(@(xy) mean(xy,1), ROIs_XY_aligned{2}, 'UniformOutput', false);
        
        REG_DATA_2_match.ROIs_XY.aligned = locfun_local_align_roi(KeyPoints.kp{1}, KeyPoints.kp{2}, centroids_2, REG_DATA_2.ROIs_XY, REG_DATA_2_match.ROIs_XY.aligned, radius_npix, min_num_refpt);
        REG_DATA_1_match.ROIs_XY.aligned = locfun_local_align_roi(KeyPoints.kp{2}, KeyPoints.kp{1}, centroids_1, REG_DATA_1.ROIs_XY, REG_DATA_1_match.ROIs_XY.aligned, radius_npix, min_num_refpt);
        
        % Get crop valid FOV: transform to binary masks and flag the ROIs that have
        % area that shrinks too much... i.e. are likely out if FOV. Remember we
        % need the ROIs that are likely common to both fields of view!
        [is_ok_bothfovs_1] = locfun_flag_outOfFov_ROIs(REG_DATA_1.ROIs_XY, REG_DATA_1_match.ROIs_XY.aligned, REG_DATA_2.imgsize, KEYPT_DATA.info_alignedImg.prepro, options.min_roi_area);
        [is_ok_bothfovs_2] = locfun_flag_outOfFov_ROIs(REG_DATA_2.ROIs_XY, REG_DATA_2_match.ROIs_XY.aligned, REG_DATA_1.imgsize, KEYPT_DATA.info_baseImg.prepro, options.min_roi_area);
        
    end


    function [REG_RES_ell] = do_match_ellipse()
        
        ROIs_XY_ell_1 = cellfun(@(roi_XY) tranform_ROI_2ellipse(roi_XY, options.ellipse_size_constrain) , REG_DATA_1.ROIs_XY, 'UniformOutput', false);
        
        if options.local_match
            
            radius_npix = 150;
            min_num_refpt = 8; % At least 3 matches are needed to provide a solution
            
            ROIs_XY_ell_2 = cellfun(@(roi_XY) tranform_ROI_2ellipse(roi_XY, options.ellipse_size_constrain) , REG_DATA_2.ROIs_XY, 'UniformOutput', false);
            
            centroids_2 = cellfun(@(xy) mean(xy,1), ROIs_XY_ell_2, 'UniformOutput', false);
            
            ROIs_XY_ell_2a = locfun_local_align_roi(KeyPoints.kp{1}, KeyPoints.kp{2}, centroids_2, ROIs_XY_ell_2, ...
                locfun_transform_roi_coord(ROIs_XY_ell_2,M,T,'2_to_1'), radius_npix, min_num_refpt);
            
        else
            ROIs_XY_ell_2a = cellfun(@(roi_XY) locfun_transform_roi_coord(...
                tranform_ROI_2ellipse(roi_XY, options.ellipse_size_constrain) ,M,T,'2_to_1') , REG_DATA_2.ROIs_XY, 'UniformOutput', false);
        end
        % Match ROIs IDs - Note that we use the image 1 as template to project all
        % ROIs, excluding the ones that are out of one of the FOVs.
        [REG_RES_ell] = matchROIs_2ses({ROIs_XY_ell_1(is_ok_bothfovs_1), ROIs_XY_ell_2a(is_ok_bothfovs_2)}, ...
            {find(is_ok_bothfovs_1), find(is_ok_bothfovs_2)}, REG_DATA_1.imgsize, options.match_opt);
    end




% -------------------------------------------------------------------------
% FINALISE RESULTS REG_RES
% -------------------------------------------------------------------------


% Add the list that was excluded
REG_RES.outofFOV_ROIs{1,1} = find(~is_ok_bothfovs_1);
REG_RES.outofFOV_ROIs{1,2} = find(~is_ok_bothfovs_2);
if any(~is_ok_bothfovs_1) || any(~is_ok_bothfovs_2)
    temp_scores = 99*ones(length(is_ok_bothfovs_1), length(is_ok_bothfovs_2));
    temp_scores(is_ok_bothfovs_1,is_ok_bothfovs_2) = REG_RES.score;
    REG_RES.score = temp_scores;
end

% Add empty lists for deleted rois
REG_RES.deleted_ROIs = {[],[]};


% Add the name of the sessions
REG_RES.sessions = {session_1, session_2};

% Unique ID of the ROIs as obtained during drawing
REG_RES.ROIsIDs = {REG_DATA_1.ROIsIDs, REG_DATA_2.ROIsIDs};


% -------------------------------------------------------------------------
% RESET ORIGINS FOR IMAGES AND COORDINATES (SYNC)
% -------------------------------------------------------------------------
% /!\  Because of transformations the coordinates might endup with negative
% values for x/y. To be able to display both images in the same coordinate
% system we have to move the origins. Hence we will need to re-adjust the
% rois and corners coordinates, as well as padding the images.

[REG_DATA_1_match, REG_DATA_2_match] = locfun_sync_origins(REG_DATA_1_match, REG_DATA_2_match);
[REG_DATA_2_match, REG_DATA_1_match] = locfun_sync_origins(REG_DATA_2_match, REG_DATA_1_match);


% -------------------------------------------------------------------------
% COMPUTE CENTROIDS
% -------------------------------------------------------------------------

% Get centroids for template data
REG_DATA_1_match.centroids.original = cellfun(@(xy) mean(xy,1), REG_DATA_1.ROIs_XY, 'UniformOutput', false);
REG_DATA_2_match.centroids.original = cellfun(@(xy) mean(xy,1), REG_DATA_2.ROIs_XY, 'UniformOutput', false);
% Get centroids for aligned data
REG_DATA_1_match.centroids.aligned = cellfun(@(xy) mean(xy,1), REG_DATA_1_match.ROIs_XY.aligned, 'UniformOutput', false);
REG_DATA_2_match.centroids.aligned = cellfun(@(xy) mean(xy,1), REG_DATA_2_match.ROIs_XY.aligned, 'UniformOutput', false);


end

function [REG_DATA_match_ref, REG_DATA_match_align] = locfun_sync_origins(REG_DATA_match_ref, REG_DATA_match_align)

% Initialise aligned sync
REG_DATA_match_align.ROIs_XY.aligned_sync = REG_DATA_match_align.ROIs_XY.aligned; 
REG_DATA_match_align.corners_XY.aligned_sync = REG_DATA_match_align.corners_XY.aligned; 

% Get the offset required for the origin of the image
addY =  floor(min(REG_DATA_match_align.corners_XY.aligned(2,[1,4]))); % [top left, bottom left, bottom right, top right]
addX =  floor(min(REG_DATA_match_align.corners_XY.aligned(1,[1,2])));

update_origin(1,2, addY)
update_origin(2,1, addX)

% Next pad the images to fill the total image
[d1_ref,d2_ref] = size(REG_DATA_match_ref.imgtiff.original_sync);
[d1_a,d2_a] = size(REG_DATA_match_align.imgtiff.aligned_sync);

if d1_ref>d1_a
    REG_DATA_match_align.imgtiff.aligned_sync = padfill_tif_img(REG_DATA_match_align.imgtiff.aligned_sync, d1_ref-d1_a, 1);
elseif d1_ref<d1_a
    REG_DATA_match_ref.imgtiff.original_sync = padfill_tif_img(REG_DATA_match_ref.imgtiff.original_sync, d1_a-d1_ref, 1);
end

if d2_ref>d2_a
    REG_DATA_match_align.imgtiff.aligned_sync = padfill_tif_img(REG_DATA_match_align.imgtiff.aligned_sync, d2_ref-d2_a, 2);
elseif d2_ref<d2_a
    REG_DATA_match_ref.imgtiff.original_sync = padfill_tif_img(REG_DATA_match_ref.imgtiff.original_sync, d2_a-d2_ref, 2);
end


    function update_origin(dim_rc,dim_xy, addZ)
        
        if addZ>0
            % Pad the aligned image
            REG_DATA_match_align.imgtiff.aligned_sync = pad_origin_tif_img(REG_DATA_match_align.imgtiff.aligned_sync, addZ, dim_rc);
            
        elseif addZ<0
            addZ = abs(addZ);
            
            % Pad the ref image
            REG_DATA_match_ref.imgtiff.original_sync = pad_origin_tif_img(REG_DATA_match_ref.imgtiff.original_sync, addZ, dim_rc);
            
            % Move the ref corners
            REG_DATA_match_ref.corners_XY.original_sync(dim_xy,:) = REG_DATA_match_ref.corners_XY.original_sync(dim_xy,:)+addZ;
            
            % Move the ref ROIs coordinates
            REG_DATA_match_ref.ROIs_XY.original_sync = cellfun(@(c) shift_coord(c, dim_xy, addZ), REG_DATA_match_ref.ROIs_XY.original_sync, 'UniformOutput', false);
            
            % Move the aligned corners
            REG_DATA_match_align.corners_XY.aligned_sync(dim_xy,:) = REG_DATA_match_align.corners_XY.aligned_sync(dim_xy,:)+addZ;
            
            % Move the aligned ROIs coordinates
            REG_DATA_match_align.ROIs_XY.aligned_sync = cellfun(@(c) shift_coord(c, dim_xy, addZ), REG_DATA_match_align.ROIs_XY.aligned_sync, 'UniformOutput', false);
            REG_DATA_match_align.ROIs_XY.globalT_sync = cellfun(@(c) shift_coord(c, dim_xy, addZ), REG_DATA_match_align.ROIs_XY.globalT_sync, 'UniformOutput', false);
            
        end
        
    end



end


function [c] = shift_coord(c, dim_xy, addZ)


c(:,dim_xy) = c(:,dim_xy)+addZ;


end

function [tpl_tiffs] = pad_origin_tif_img(tpl_tiffs, z_add, dim)

if dim==2
    % add columns at left
    tpl_tiffs = cat(2,zeros(size(tpl_tiffs,1),z_add),tpl_tiffs);
else
    % add rows at top
    tpl_tiffs = cat(1,zeros(z_add, size(tpl_tiffs,2)),tpl_tiffs);
end

end

function [tpl_tiffs] = padfill_tif_img(tpl_tiffs, z_add, dim)

if dim==2
    % add columns at right
    tpl_tiffs = cat(2,tpl_tiffs,zeros(size(tpl_tiffs,1),z_add));
else
    % add rows at bottom
    tpl_tiffs = cat(1,tpl_tiffs,zeros(z_add, size(tpl_tiffs,2)));
end

end

% function [] = locfun_adjust_matchedROIs(matched_ROIs)
% % Get correct alignment by find the best correlation patch using the
% % aligned data
% patch_halfsiz = 20;
% patch_xy_shift = 0;% 0 # pixels
%
% for iMatch = 1:size(matched_ROIs,2)
%
%     % Get image patch

%     % Do correlation - but just translate image
%
%     % If good correlation, run the rest - if not don't check it
%
% end
%
% end




function [ROIs_XY_aligned] = locfun_local_align_roi(kp1, kp2, centroids_2, ROIs_XY_2, ROIs_XY_aligned, radius_npix, min_num_refpt)

num_kp = size(kp2,1);% [nPts x 2]

for iRoi=1:length(ROIs_XY_2)
    
    roi_c = centroids_2{iRoi};
    
    % Find closest KP 
    
    dist_roi_kp = sqrt(sum((repmat(roi_c, [num_kp, 1]) - kp2).^2, 2));
    is_kp_close = dist_roi_kp < radius_npix;
    
    % If not enough ref points, do nothing
    %     if (sum(is_roi_close)+sum(is_kp_close)) < min_num_refpt
    if (sum(is_kp_close)) < min_num_refpt
        continue
    end

    % Get affine transformation using local KP and matched rois
    [M, T] = locfun_get_transformation_KP(kp1(is_kp_close,:),kp2(is_kp_close,:));
    % size T: [2 x 1]
    
    % Rois coordinates are [Npts x 2(x,y)]
    
    % Transfo ROIs coordinates from session 2 to coord system of session 1
    
    roi_xy_alg = M*ROIs_XY_2{iRoi}'+ repmat(T, [1, size(ROIs_XY_2{iRoi},1)]);
    
    ROIs_XY_aligned{iRoi} = roi_xy_alg';
    
end
end


function [is_ok_bothfovs] = locfun_flag_outOfFov_ROIs(ROIs_XY, ROIs_XY_aligned, imgsize, img_prepro, min_roi_area)
% d1 = rows = y
% d2 = columns = x

d1 = imgsize(1);
d2 = imgsize(2);


% Check boundaries (default: min=1 max=d, but might have cropped image)

if img_prepro.do_crop
    x_min=double(img_prepro.crop_edges.columns(1)); y_min=double(img_prepro.crop_edges.rows(1));
    x_max=d2-double(img_prepro.crop_edges.columns(2)); y_max=d1-double(img_prepro.crop_edges.rows(2));
else
    x_min=1; y_min=1;
    x_max=d2; y_max=d1;
end

% Find ROIs that have coordinates out of bounds
has_coord_out_alg = cellfun(@(R)any(R(:,1)<x_min | R(:,1)>x_max) | any(R(:,2)<y_min | R(:,2)>y_max) , ROIs_XY_aligned);
has_coord_out = has_coord_out_alg | cellfun(@(R)any(R(:,1)<x_min | R(:,1)>x_max) | any(R(:,2)<y_min | R(:,2)>y_max) , ROIs_XY);


if any(has_coord_out)
    % Original area for each ROI
    A_orig = cellfun(@(R)polyarea(R(:,1),R(:,2)), ROIs_XY(has_coord_out), 'UniformOutput', true);
    
    % A_orig_min = min(A_orig);
    
    % Area after cropped image
    A_orig_crop = cellfun(@(R)polyarea( ...
        min(max(R(:,1),x_min),x_max), min(max(R(:,2),y_min),y_max) ), ROIs_XY(has_coord_out), 'UniformOutput', true);
    
    % New area for each ROI, after projection to other FOV - Bound aligned
    % coordinates to min=1 max=d
    A_new = cellfun(@(R)polyarea( ...
        min(max(R(:,1),x_min),x_max), min(max(R(:,2),y_min),y_max) ), ROIs_XY_aligned(has_coord_out), 'UniformOutput', true);
    
    % Falg the areas =0 or < min_roi_area *A_orig
    
    is_out_of_fov = false(size(ROIs_XY));
    is_out_of_fov(has_coord_out) = (A_orig_crop < min_roi_area *A_orig) | A_new < min_roi_area *A_orig ;
    % is_out_of_fov(has_coord_out) = A_new < min_roi_area *A_orig_min;
    
    is_ok_bothfovs = ~is_out_of_fov;
    
else
    is_ok_bothfovs = true(size(ROIs_XY));
end


end




function [Cxy, Cxy_aligned] = locfun_transform_imgcorners(M, T, Cxy, flag)
% d1 d2 = size of image (row/column)
% T = [tx, ty]'
% NOT confusing at all, nor prone to bugs...

% image: [top left, bottom left, bottom right, top right]

% top left image, aka origin (x=0.5,y=0.5)
% bottom left image, aka (x=0.5,y=d1+0.5)
% bottom right image, aka (x=d2+0.5,y=d1+0.5)
% top right image, aka (x=d2+0.5,y=0.5)

% d1 = imgsize(1);
% d2 = imgsize(2);
% Cxy = [[0.5;0.5] , [0.5;d1+0.5], [d2+0.5;d1+0.5], [d2+0.5;0.5]];

Cxy_aligned = zeros(2,4);

switch flag
    case '2_to_1'
        
        Cxy_aligned = M*Cxy + T;
        
        
    case '1_to_2'
        
        Cxy_aligned = M \ (Cxy - T);
        
end


end


function [ROIs_XY_aligned] = locfun_transform_roi_coord(ROIs_XY,M,T,flag)

% KP1 = (M * KP2' + T)'   ->  KP1 = KP2 * M'+ T'     % [nPts x 2]
% KP2 = (M \ (KP1' - T))'   ->  KP2 =  (KP1 - T') * inv(M')

% Transpose
M = M';
T = T';

switch flag
    case '2_to_1'
        
        % Transfo ROIs coordinates from session 2 to coord system of
        % session 1
        ROIs_XY_aligned = cellfun(@(roi) roi*M + repmat(T, [size(roi,1),1]), ROIs_XY, 'UniformOutput', false);
        
        
    case '1_to_2'
        
        % Transfo ROIs coordinates from session 1 to coord system of
        % session 2
        ROIs_XY_aligned = cellfun(@(roi) (roi-repmat(T, [size(roi,1),1]))/M, ROIs_XY, 'UniformOutput', false);
        
end


end


function [M, T] = locfun_get_transformation(kp1,kp2,rois_c_1,rois_c_2)

% keypoints coordinates
% kp1, kp2 [nPts x 2]

refpt_1 = [kp1;rois_c_1]; refpt_2 = [kp2;rois_c_2];
% Get the affine parameters (we transform Kp2 back to image 1)
A = locfun_get_transfoMat(refpt_2);

% Solve A*z=b -> z = (A'A)^-1 * A' * b

b=reshape(refpt_1',[2*size(refpt_1,1),1]);

[M, T] = locfun_get_affineparam(A, b);

end


function [M, T] = locfun_get_transformation_KP(kp1, kp2)

% Get the affine parameters (we transform Kp2 back to image 1)
A = locfun_get_transfoMat(kp2);

% Solve A*z=b -> z = (A'A)^-1 * A' * b

b=reshape(kp1',[2*size(kp1,1),1]);

[M, T] = locfun_get_affineparam(A, b);

end

function [M, T] = locfun_get_affineparam(A, b)

z=pinv(A)*b;

% Extract parameters

M = [z(1),z(2) ; z(3),z(4)];
T = [z(5);z(6)];

end


function [A] = locfun_get_transfoMat(kp2)

% Affine transformation

Ac = arrayfun(@(i)[kp2(i,:), [0,0,1,0] ; [0,0], kp2(i,:), [0,1]], 1:size(kp2,1), 'UniformOutput', false);
A=cat(1,Ac{:});  % [nPts*2 x 6]

% Solve :  KP1 = M * KP2 + T
% [ kp1_x ,  = [ m1 m2 ,  [ kp2_x ,  + [ tx ,
%   kp1_y ]      m3 m4 ]    kp2_y ]      ty ]
%
% Tranform to : A * z = KP1
% where A =
% [ x y 0 0 1 0 ,
%   0 0 x y 0 1 ,
%   ...          ]
% (x,y are the coord of kp2)
% and Z = [ m1 m2 m3 m4 tx ty ]'

end


function [REG_RES] = locfun_merge_REGRES(REG_RES, REG_RES_ell)

% COMBINE THE RESULTS OF REGISTRATION

% Find the pairs that already exist
ellipse_pair_exists = arrayfun(@(n1,n2) any((n1==REG_RES.matched_pairs(1,:))&(n2==REG_RES.matched_pairs(2,:))), ...
    REG_RES_ell.matched_pairs(1,:), REG_RES_ell.matched_pairs(2,:));

% Find the new pairs registered in ELLIPSE, without conflict with existing
% pairs (i.e. both rois must be in the nonmatched group):
ellipse_pair_is_valid = arrayfun(@(n1,n2, is_new) is_new && (any(n1==REG_RES.nonmatched_ROIs{1}) && any(n2==REG_RES.nonmatched_ROIs{2})), ...
    REG_RES_ell.matched_pairs(1,:), REG_RES_ell.matched_pairs(2,:), ~ellipse_pair_exists);


% Add them to the main REG_RES
REG_RES.matched_pairs = cat(2, ...
    REG_RES.matched_pairs, ...
    REG_RES_ell.matched_pairs(:,ellipse_pair_is_valid));

% Remove corresponding ROIs from nonmatched
rmv_nonmatch_1 = arrayfun(@(n1) any(n1==REG_RES_ell.matched_pairs(1, ellipse_pair_is_valid)), REG_RES.nonmatched_ROIs{1});
rmv_nonmatch_2 = arrayfun(@(n2) any(n2==REG_RES_ell.matched_pairs(2, ellipse_pair_is_valid)), REG_RES.nonmatched_ROIs{2});
REG_RES.nonmatched_ROIs{1} = REG_RES.nonmatched_ROIs{1}(~rmv_nonmatch_1);
REG_RES.nonmatched_ROIs{2} = REG_RES.nonmatched_ROIs{2}(~rmv_nonmatch_2);


% Possible other outputs:
% M_union = [M{1},M{2}(:,HMatch_results.nonmatched_2)];
% % % R = sparse(matched_ROIs(:,1),matched_ROIs(:,2),1,K1,K2);


end

