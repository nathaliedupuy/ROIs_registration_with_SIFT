% INFO ABOUT DATA STRUCTURES
%
% * KeyPoints data structure: Generated after run alignement with SIFT
%       filename format = KEYPT_DATA_[session_2align]_to_[session_template].mat
%       structure with fields:
%           - params_match
%           - info_alignedImg : (see below)
%           - info_baseImg : (see below)
%           - num_goodcorr_matches : number of matches found with good
%           correlation, before rmv outliers using MahalD
%           - matches : cell array with all the feature matches found
%
%       'info_alignedImg' and 'info_baseImg' are structures with fields:
%
%           - prepro : preprocessing info - struct with fields:
%               * shift_xy = [+x, +y] / how much was cropped left (columns)
%               and top (rows) (used to shift coordinates space after
%               feature matching)
%               * do_crop : boolean / if template image was cropped for
%               feature matching
%               * crop_edges (optional, if do_crop=True)
%
%           - name : name of image template that was used
%           - SIFT_kpts : total number of keypoints found with SIFT
%           - size : original image size
%
%       'matches' is a structure with fields:
%           - MahalD : Mahalanobis distance
%           - pcorr : pearson correlation score
%           - kp1_xy = [x, y, orientation] of keypoint descriptor for
%           session_template
%           - kp2_xy = [x, y, orientation] of keypoint descriptor for
%           session_2align
%           - is_edg_patch : boolean - whether the kp was close to the
%           edge of the image
%
%
%
% * REG_RES is the structure with the resgistation results.
%
%   - sessions = {sesname1 , sesname2}
%   - score         [n1 x n2]
%           distance between rois
%           n1=nb rois in sesname1 and n2=nb rois in sesname2
%   - matched_pairs  [2 x m]
%           indices of the m pairs of rois that were matched without
%           conflicts across sessions
%   - nonmatched_ROIs    {[1 x nm1] , [1 x nm2]}
%           indices of the rois not paired
%   - conflicting_chain [2 x ccm]
%           pairs of matched rois which are not transitive across sessions
%           (e.g. match: roiA-roiB, roiA-roiC, but not roiB-roiC)
%   - incomplete_chain [2 x icm]
%           pairs of matched rois which are not transitive across sessions,
%           but where there is no conflict (e.g. match: roiA-roiB,
%           roiA-roiC, but roiB and roiC are both in nonmatched)
%   - missed_pairs: [2 x mp]
%           pairs of rois that could be a potential match as part of an
%           incomplete chain (see above)
%   - outofFOV_ROIs      {[1 x nof1] ,  [1 x nof2]}
%           indices of the rois that are, at least in part, out of one or
%           more fields of view
%   - deleted_ROIs: {[]  []}
%           empty for now; this field can be edited later with the GUI
%   - ROIsIDs       {[1 x n1]  [1 x n2]}
%           real roi IDs, the ones created in imageJ
%
% /!\ The rois indices in the lists 'matched_pairs', 'nonmatched' and
% 'outofFOV' are not the roi IDs. You can access their IDs in ROIsIDs(i).

function [REG_RES, REG_DATA, KeyPoints] = run_ROIsRegistration(path_regdata, registration_UID, list_sessions, options)

path_regdata = fullfile(path_regdata, registration_UID);

num_ses = length(list_sessions);

% First prep images for each session
if options.prep_images
    
    fprintf('Prepare sessions template images for registration %s \n', registration_UID)
    if options.do_enhance_template
        the_inputs.path_source = path_regdata;
        fprintf('\tTemplate images will be enhanced  \n')
        fprintf('\tsession')
    end
    for iSes = 1:num_ses
        
        if options.do_enhance_template
            
            if iSes==num_ses
                fprintf('_%i/%i\n', iSes, num_ses)
            else
                fprintf('_%i/%i', iSes, num_ses)
            end
            
            the_inputs.name_source = [list_sessions{iSes},'_original'];
            the_inputs.name_output = list_sessions{iSes};
            
            % Run image enhance, and Scale to [0,255] and convert to uint8
            cmd = get_python_cmd('enhance_template', options.root_path, options.python_path, options.settingsfile, the_inputs);
            system(cmd, '-echo');
            
        else
            
            % Open original template image
            TifLink = Tiff(fullfile(path_regdata, [list_sessions{iSes},'_original.tif']), 'r');
            template_img = TifLink.read();
            TifLink.close()
            % Scale to [0,255] and convert to uint8
            template_tiff_file = fullfile(path_regdata,[list_sessions{iSes},'.tif']);
            save_tifftemplate(template_tiff_file, template_img)
            
        end
        
    end
    disp('...Prep images done.')
end

% Align the field of views using SIFT (feature detection and matching)
if options.align_fovs
    
    the_inputs = struct();
    the_inputs.path_regdata = path_regdata;
    
    fprintf('Align sessions template images for registration %s \n', registration_UID)
    
    
    for iSes = 1:num_ses-1
        
        % Choose a base template and align the other sessions to it
        ses2align = iSes+1:num_ses;
        
        the_inputs.expID_base = list_sessions{iSes};
        the_inputs.expIDs_2align = list_sessions(ses2align);
        
        % Run feature extraction and matching with SIFT:
        cmd = get_python_cmd('match_kpts_SIFT', options.root_path, options.python_path, options.settingsfile, the_inputs);
        system(cmd, '-echo');
        
    end
    
    
    
    disp('...Image alignment done.')
    
end

% Register the ROIs

if options.register_rois
    
    fprintf('Register rois for project %s \n', registration_UID)
    
    % For each session, load template image, rois coordinates, etcc..
    REG_DATA = struct(...
        'session', list_sessions(:), ...
        'data', cell(num_ses,1), ...
        'aligned_data' , cell(num_ses,1));
    %     % Init blank structure for aligned data
    %     banana_struct = [list_sessions(:)';cell(1,num_ses)];
    for iSes = 1:num_ses
        
        % Load template image
        TifLink = Tiff(fullfile(path_regdata, [list_sessions{iSes},'.tif']), 'r');
        REG_DATA(iSes).data.imgtiff = TifLink.read();
        TifLink.close()
        % Save size of image
        REG_DATA(iSes).data.imgsize = size(REG_DATA(iSes).data.imgtiff);
        % Image Corners
        d1 = REG_DATA(iSes).data.imgsize(1);
        d2 = REG_DATA(iSes).data.imgsize(2);
        REG_DATA(iSes).data.corners_XY = [[0.5;0.5] , [0.5;d1+0.5], [d2+0.5;d1+0.5], [d2+0.5;0.5]];
        % Load ROIs coordinates
        rois_coord = load(fullfile(path_regdata, [list_sessions{iSes},'_rois.mat']));
        REG_DATA(iSes).data.ROIs_XY = rois_coord.ROIs_XY;
        REG_DATA(iSes).data.ROIsIDs = rois_coord.ROIsIDs;
        % Init aligned data
        REG_DATA(iSes).aligned_data = struct();
        %         ses_banana_struct = banana_struct(:,~ismember(list_sessions,list_sessions{iSes}));
        %         REG_DATA(iSes).aligned_data = struct(ses_banana_struct{:});
        
    end
    
    
    % Note for REG_DATA(iSes).aligned_data.(sesname): We have to add a
    % character at the start of session name to save as field in the
    % structure, in case the session name starts with a number.
    
    num_pairs = num_ses*(num_ses-1)/2;
    REG_RES = cell(num_pairs,1);
    KeyPoints = cell(num_pairs,1);
    
    counter_pair=0;
    for iSes = 1:num_ses-1
        
        session_1 = list_sessions{iSes};
        
        % Run registraion between 2 sessions
        for iSes2=iSes+1:num_ses
            counter_pair=counter_pair+1;
            
            session_2 = list_sessions{iSes2};
            
            fprintf('\tRegister rois %s and %s \n', session_1, session_2)
            
            [REG_RES{counter_pair}, ...
                REG_DATA(iSes).aligned_data.(['s',session_2]), REG_DATA(iSes2).aligned_data.(['s',session_1]), ...
                KeyPoints{counter_pair}] = registerROIs_2sessions(path_regdata, session_1, REG_DATA(iSes).data, session_2, REG_DATA(iSes2).data, options.roireg);
        end
        
    end
    
    REG_RES = cat(1,REG_RES{:});
    
    if num_ses>2
        fprintf('\tChecking matches are consistent across sessions \n')
        % Check that there are no conflicting matches
        REG_RES = check_reg_across_sessions(list_sessions, REG_RES);
        
        % Quick we did not add/loose rois on the way...
        [no_roi_missing] = check_no_roi_is_missing(REG_RES);
        if ~no_roi_missing
            error('Rois got lost on the way during check registration across sessions')
        end
    end
    
    
    
    disp('...ROIs Registration done.')
    
else
    
    REG_RES=[];
    REG_DATA=[];
    KeyPoints=[];
    
end




end


function save_tifftemplate(template_tiff_file, template_img)

[d1, d2] = size(template_img);

tifftagstruct = struct(...
    'Compression',Tiff.Compression.None,...
    'ImageWidth',d2,'ImageLength',d1,...
    'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky, ...
    'BitsPerSample', 8, ... % 'uint8', 'int8' (choose with 'SampleFormat')
    'SampleFormat', Tiff.SampleFormat.UInt, ...
    'SamplesPerPixel', 1, ... % for greyscale
    'Photometric',Tiff.Photometric.MinIsBlack, ...
    'Software','MATLAB');

t = Tiff(template_tiff_file, 'w');
t.setTag(tifftagstruct);
t.write(uint8(round(minmaxscale255(template_img))));
t.close();

end

function [img] = minmaxscale255(img)


img_min = min(img(:));

img_max = max(img(:));

img =(img - img_min) / (img_max - img_min);

img = img*255;

end

function [no_roi_missing] = check_no_roi_is_missing(REG_RES)

no_roi_missing = arrayfun(@(R)f_check_siz(R), REG_RES);

no_roi_missing = all(no_roi_missing);

end

function [no_roi_missing] = f_check_siz(REG_RES_n)

no_roi_missing = true;
for r=1:2
    if no_roi_missing
        no_roi_missing = size(REG_RES_n.matched_pairs,2) +...
            size(REG_RES_n.nonmatched_ROIs{r},2) +...
            size(REG_RES_n.outofFOV_ROIs{r},2) + ...
            size(REG_RES_n.deleted_ROIs{r},2) + ...
            size(REG_RES_n.conflicting_chain,2) + ...
            size(REG_RES_n.incomplete_chain,2) + ...
            size(REG_RES_n.missed_pairs,2) == length(REG_RES_n.ROIsIDs{r});
    end
end

end
