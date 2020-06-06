function prepare_data_nrlab(path_regdata, registration_UID, path_workingdir, the_animal, the_fov, list_sessions, list_sessions_newnames, options)
% The template is the first recording by default, and hence the ROIs mask
% taken is the one corresponding to the first recording.

path_regdata = fullfile(path_regdata, registration_UID);

if ~iscell(list_sessions)
    list_sessions={list_sessions};
end

num_ses = length(list_sessions);

% Check path to registation data exists
if ~exist(path_regdata, 'dir')
    mkdir(path_regdata)
    fprintf('Create directory %s \n', path_regdata)
end

% Check options for motion correction
if options.do_align_across &&  ~isfield(options, 'mc_max_shift')
    options.mc_max_shift = 50;
end


% For each session, load desired recordings and align all trials (rigid MC).
fprintf('Prepare data for registration %s \n', registration_UID)
if options.print_mc_msg
    msg2print='session_%i/%i \n';
else
    msg2print='_%i/%i';
    fprintf('session')
end
for iSes = 1:num_ses
    
    if isempty(the_fov)
        the_roiset_folder = dir(fullfile(path_workingdir, list_sessions{iSes}, 'ROIs', ['ROISets_',the_animal,'_*']));
        original_tifffolder = fullfile(path_workingdir, list_sessions{iSes}, 'ROIs', the_roiset_folder.name);
    else
        original_tifffolder = fullfile(path_workingdir, list_sessions{iSes}, 'ROIs', ['ROISets_',the_animal,'_',the_fov]);
    end
    
    
    if ~options.print_mc_msg && (iSes==num_ses)
        fprintf([msg2print,'\n'], iSes, num_ses)
    else
        fprintf(msg2print, iSes, num_ses)
    end
    
    
    % Prep the recordings that will be included in the template
    ses_recs_tif = get_all_recs_tiff(original_tifffolder);
    num_recs = length(ses_recs_tif);
    
    if options.export_images
        % Open all the average tiff images of all recordings.
        % Get dimensions using the first recording
        tiffInfo = imfinfo(fullfile(original_tifffolder, ses_recs_tif{1}));
        d2 = tiffInfo(1).Width;
        d1 = tiffInfo(1).Height;
        Y_recs = zeros(d1,d2,num_recs,'single');
        for iRec=1:num_recs % NO NEED TO PARFOR IT'S SUPER FAST
            
            TifLink = Tiff(fullfile(original_tifffolder, [ses_recs_tif{iRec}]), 'r');
            
            Y_recs(:,:,iRec) = TifLink.read();
            
            TifLink.close()
            
        end
        
        % Align to first rec by default, or simply average if the tiffs
        % were already aligned across trials. Note that the alignement will
        % align to the first tiff image (template by default)
        if options.do_align_across
            
            %         fprintf('\t Align rec images with normcorre (Pnevmatikakis, 2016) \n')
            
            options_rigid = NoRMCorreSetParms('d1',d1,'d2',d2,'max_shift',options.mc_max_shift,'us_fac',50,'upd_template',false, ...
                'print_msg', options.print_mc_msg);
            
            % We a custom lighter version of normcorre
            [M1] = normcorre_rigid_light(Y_recs(:,:,2:end), options_rigid, Y_recs(:,:,1));
            
            template_img = mean(cat(3, Y_recs(:,:,1), M1), 3);
            
            
        else
            
            if num_recs>1
                template_img = mean(Y_recs, 3);
            else
                template_img = Y_recs;
            end
            
        end
        
        
        % Save original template before any modification
        template_tiff_file = fullfile(path_regdata,[list_sessions{iSes},'_original.tif']);
        save_tifftemplate_original(template_tiff_file, template_img, d1, d2)
        
    end
    
    % Copy the ROIzip file
    if options.export_rois
        
        [~,name_firstrec,~] = fileparts(ses_recs_tif{1});
        first_rec_rois = [name_firstrec,'.zip'];
        
        % Get the ROIs coordinates (ImageJ -> MATLAB function ReadImageJROI by Dylan Muir)
        [ROIs_XY, ROIsIDs] = load_ROIs_XY(fullfile(original_tifffolder, first_rec_rois), options.rois_2exclude);
        % size: [1xN]
        
        % Save coordinates and IDs
        ses_rois_coord = [list_sessions{iSes},'_rois.mat'];
        save(fullfile(path_regdata, ses_rois_coord), 'ROIs_XY', 'ROIsIDs')
        
    end
    
    
end
disp('...Prep data done.')

end

function save_tifftemplate_original(template_tiff_file, template_img, d1, d2)

tifftagstruct = struct(...
    'Compression',Tiff.Compression.None,...
    'ImageWidth',d2,'ImageLength',d1,...
    'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky, ...
    'BitsPerSample', 32, ... % 'uint32', 'int32'
    'SampleFormat', Tiff.SampleFormat.IEEEFP, ...
    'SamplesPerPixel', 1, ... % for greyscale
    'Photometric',Tiff.Photometric.MinIsBlack, ...
    'Software','MATLAB');

t = Tiff(template_tiff_file, 'w');
t.setTag(tifftagstruct);
t.write(template_img);
t.close();

end


function [list_sesRecs] = get_all_recs_tiff(path_tiffs)

% Find all tiffs in the directory
listing_tif = dir(fullfile(path_tiffs, '*.tif'));
listing_tif2 = dir(fullfile(path_tiffs, '*.tiff'));

list_sesRecs={};
if ~isempty(listing_tif)
    list_sesRecs = [list_sesRecs, {listing_tif.name}];
end
if ~isempty(listing_tif2)
    list_sesRecs = [list_sesRecs, {listing_tif2.name}];
end

% listing = dir(path_tiffs);
% recs = {listing(~[listing.isdir]).name}; is_tiff = cellfun(@(name)strcmp(name(end-3:end),'.tif'), recs);
% list_SesRecs = recs(is_tiff);


end


function [ROIs_XY, cellsIDs] = load_ROIs_XY(roiszipfile, rois_2exclude)
% Loads the ROI spatial coordinates

% ATTENTION HAS TO BE strType: 'Polygon'

sROI = ReadImageJROI(roiszipfile);
% ReadImageJROI (Dylan Muir)
% https://github.com/DylanMuir/ReadImageJROI/blob/master/ReadImageJROI.m

cellsIDs = 1:length(sROI);
% Note that we could extract their real ID (see last line of code)

if ~isempty(rois_2exclude)
    idx_rois_2keep = setdiff(1:length(sROI), rois_2exclude);
    sROI = sROI(idx_rois_2keep);
    cellsIDs = cellsIDs(idx_rois_2keep);
end

% Also remove if other than Polygon...
is_polygon = cellfun(@(s) strcmp(s.strType, 'Polygon'), sROI, 'UniformOutput', true);
if any(~is_polygon)
    %     fprintf('\n')
    warning(['More non-polygon ROIs found... removing them (n=',num2str(sum(~is_polygon)),')'])
    sROI = sROI(is_polygon);
    cellsIDs = cellsIDs(is_polygon);
end

ROIs_XY = cellfun(@(s) s.mnCoordinates, sROI, 'UniformOutput', false);

% Note that we could extract their real ID:
% cellsIDs = cellfun(@(s) s.strName, sROI, 'UniformOutput', false);

end

