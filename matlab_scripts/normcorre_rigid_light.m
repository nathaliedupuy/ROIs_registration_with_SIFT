function [M_final,shifts,template,options,col_shift] = normcorre_rigid_light(Y,options,template)

% Light version of 'normcorre' (nathalie dupuy)
%   See normcorre for original + authors
%       https://github.com/flatironinstitute/NoRMCorre
%
% **EDIT** 
%   - skip filetype check, as we know it is tiff
%   - skip dimension check, as we know it is not a volume, but 2D image
%   - skip set default parameters if not present - we have to give it
%   - skip plot_flag
%   - remove grid: we do not use patches (rigid only here)
%   - remove mat2cell_ov as it is also obsolete here (no grid)
%   - remove FFT shift method option
%   - remove add value
%   - IMPORTANT: changed fill value with median pixel value, unless user requests NaN

% ***** Clarification about bin/buffers parameters:
% options.tiff.mc.normcorr.init_batch = 200; % length of initial batch to get template (default: 100)(their demo: 200) [number of frames]
% options.tiff.mc.normcorr.upd_template = true; % flag for online template updating (default: true)
% %   Next ones are only relevant if 'upd_template' is true:
% %   >   Use a buffer to estimate a new temporary template, and use it to
% %       update the current template.
% options.tiff.mc.normcorr.bin_width = 200; % width of each bin for buffer (default: 200) [number of frames]
% %   >   Select method for tpl update:  
% %           method{1}: avg/median previous tpl estimate with the new one
% %           method{2}: avg/median the buffer to get the new temp template
% options.tiff.mc.normcorr.method= {'median';'mean'}; % method for averaging the template (default: {'median';'mean'})
% %   >   If method{1} is 'median' :
% options.tiff.mc.normcorr.buffer_width = 50; % number of templates to keep in memory (default: 50) [number of frames]



% -------------------------------------------------------------------------
% INPUTS
% Y:                Input data, can be already loaded in memory as a 3D
%                   tensor, a memory mapped file, or a pointer to a tiff stack
% options:          options structure for motion correction (optional, rigid registration is performed if not provided)
% template:         provide template (optional)

% OUTPUTS
% M_final:          motion corrected data
% shifts:           originally calculated shifts
% template:         calculated template
% options:          options structure (if modified)
% col_shift:        relative shift due to bi-directional scanning



% Check if did single precision
if ~isa(Y,'single')
    Y = single(Y);
end
% Dimension is 2 (2D image x nframes)
[d1,d2,num_frames] = size(Y);


% Overwrite method to CUBIC
options.shifts_method = 'cubic';
    

%% Check options are valid

options.bin_width = min(options.bin_width,num_frames+1);
options.mem_batch_size = min(options.mem_batch_size,num_frames);
options.buffer_width = min(options.buffer_width,ceil(num_frames/options.bin_width));
options.init_batch = min(options.init_batch,num_frames);

init_batch = options.init_batch;
tpl_updt_method = options.method;
max_shift = options.max_shift;

    

%% Check for offset due to bi-directional scanning

if isempty(options.col_shift)
    if options.correct_bidir
        col_shift = correct_bidirectional_offset(Y,options.nFrames,options.bidir_us);
    else
        col_shift = 0;
    end
else
    col_shift = options.col_shift;
end
options.col_shift = col_shift;
if col_shift~=0
    if options.print_msg; fprintf('Offset %1.1d pixels due to bidirectional scanning detected. \n',col_shift); end
%     if strcmpi(options.shifts_method,'fft')
%         options.shifts_method = 'cubic';
%         if options.print_msg; fprintf('Cubic shifts will be applied. \n'); end
%     end
end

%% Read initial batch and compute template

if nargin < 3 || isempty(template)
    init_batch = min(num_frames,init_batch);
    interval = ceil(num_frames/2-init_batch/2+1):floor(num_frames/2+init_batch/2);
    if options.print_msg; fprintf('Registering the first %i frames just to obtain a good template....',init_batch); end
    Y_temp = Y(:,:,interval);
    template = median(Y_temp,3);
    fftTemp = fftn(template);
    for t = 1:size(Y_temp,3)
        [~,Greg] = dftregistration_min_max(fftTemp,fftn(Y_temp(:,:,t)),options.us_fac,-max_shift,max_shift,options.phase_flag);
        M_temp = real(ifftn(Greg));
        template = template*(t-1)/t + M_temp/t;
    end
    if options.print_msg; fprintf('..done. \n'); end
else
    if ~isa(template,'single')
        template = single(template);
    end
end
fftTemp = fftn(template);


%% Prep buffer and output

if options.upd_template
    buffer = nan(d1,d2,options.bin_width,'single');
    if strcmpi(tpl_updt_method{1},'median')
        buffer_med = nan(d1,d2,options.buffer_width,'double'); 
        % Note: it was not initialised in the original code; when created
        % it was then a double. This is important if wish to replicate the
        % results of the original code.
    end
end

M_final = zeros([d1,d2,num_frames],'single'); 
shifts = zeros(1,1,2,num_frames);


if options.print_msg; fprintf('Template initialization complete.  Now registering all the frames with new template. \n'); end


%% Start MC loop over frames

if strcmpi(options.boundary,'nan')
    fill_value = NaN;
else
    fill_value = median(Y(:));
end

cnt_buf = 0;
prevstr = [];
lb = -max_shift;
ub = max_shift;

for it = 1:options.iter
    
    shifts_temp = zeros(1,1,2);
    
    for t = 1:num_frames
        
        Yt = Y(:,:,t);
        
        minY = min(Yt(:));
        maxY = max(Yt(:));

        fftY = fftn(Yt);
      
        [output,Greg] = dftregistration_min_max(fftTemp,fftY,options.us_fac,lb,ub,options.phase_flag);
        
        M_temp = real(ifftn(Greg));
        M_temp = remove_boundaries(M_temp,output(3:end),'copy',template);
        
        if options.upd_template
            ind = rem(t,options.bin_width) + options.bin_width*(rem(t,options.bin_width)==0);
            buffer(:,:,ind) = M_temp;
        end
        
        shifts_temp(1) = output(3);
        shifts_temp(2) = output(4);
        
        shifts(1,1,:,t) = shifts_temp;        
   
        % NOTE: removed options.shifts_method : 'fft'
        
        % NOTE: both do not give the EXACT numerical same
        % results... but technically should be the same. To be
        % conservative, let's use the original code.
        shifts_up = imresize(shifts_temp,[d1,d2]);
        shifts_up(2:2:end,:,2) = shifts_up(2:2:end,:,2) + col_shift;
        Mf = imwarp(Yt,-cat(3,shifts_up(:,:,2),shifts_up(:,:,1)),options.shifts_method,'FillValues',fill_value);
        %                 shifts_up = repmat(shifts_temp,[d1,d2]);
        %                 shifts_up(2:2:end,:,2) = shifts_up(2:2:end,:,2) + col_shift;
        %                 Mf = imwarp(Yt,-cat(3,shifts_up(:,:,1,2),shifts_up(:,:,1,1)),options.shifts_method,'FillValues',fill_value);
        Mf(Mf<minY) = minY;
        Mf(Mf>maxY) = maxY;

        M_final(:,:,t) = Mf;
        

        if options.upd_template && mod(t,options.bin_width) == 0 
            if options.print_msg
                str=[num2str(t), ' out of ', num2str(num_frames), ' frames registered, iteration ', num2str(it), ' out of ', num2str(options.iter), '..'];
                refreshdisp(str, prevstr, t);
                prevstr=str; 
                %fprintf('%i out of %i frames registered, iteration %i out of %i \n',t,T,it,options.iter)
            end
            % ** EDIT: get new slice template from avg/median buffer (options.bin_width)       
            if strcmpi(tpl_updt_method{2},'mean')
                new_temp = nanmean(buffer,3);
            elseif strcmpi(tpl_updt_method{2},'median')
                new_temp = nanmedian(buffer,3);
            end
            
            % ** EDIT: update template from new slice template (options.bin_width)
            % + previous buffer template (options.buffer_width)
            if strcmpi(tpl_updt_method{1},'mean')
                cnt = t/options.bin_width + 1;
                template = template*(cnt-1)/cnt + new_temp*1/cnt;
                
            elseif strcmpi(tpl_updt_method{1},'median')
                cnt_buf = cnt_buf + 1;
                if cnt_buf <= options.buffer_width
                    buffer_med(:,:,cnt_buf) = new_temp;
                else
                    buffer_med = circshift(buffer_med,[zeros(1,2),-1]);
                    buffer_med(:,:,options.buffer_width) = new_temp;
                end
                template = nanmedian(buffer_med,3);
            end
            fftTemp = fftn(template);
            
        end        

    end

if options.print_msg; fprintf('\n'); end
    
if options.print_msg; fprintf('done. \n'); end

end