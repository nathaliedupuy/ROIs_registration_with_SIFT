function [cmd] = get_python_cmd(cmd_name, root_path, python_path, settingsfilename, the_inputs)

% Important! python command call breaks with whitespace - it separates
% strings, so implemented correction of whitespaces:
[root_path, python_path] = check_whitespaces_in_path(root_path, python_path);

% if iscell(paths_data)
%     paths_data = cellfun(@(p)replace(p,' ','\ '), paths_data, 'UniformOutput', false);
% else
%     paths_data = replace(paths_data,' ','\ ');
% end
  

switch cmd_name
        
    case 'enhance_template'
        
        path_source = replace(the_inputs.path_source,' ','\ ');
        
        cmd = sprintf('%s %s %s %s %s "" %s %s', python_path, fullfile(root_path,'python_scripts','enhance_fov_img.py'), ...
            fullfile(root_path, settingsfilename), ...
            path_source, the_inputs.name_source, the_inputs.name_output, '1');
        
    case 'match_kpts_SIFT'
        
        path_regdata = replace(the_inputs.path_regdata,' ','\ ');
        
        num_exp2align = length(the_inputs.expIDs_2align);
        list_s = repmat(' %s', 1, num_exp2align);
        
%         cmd = sprintf('%s %s', python_path, fullfile(root_path,'test_install','dummy1.py'));
        
        cmd = sprintf(['%s %s %s %s %s', list_s], python_path, fullfile(root_path,'python_scripts','match_fov_kpts_SIFT.py'), ...
            fullfile(root_path, settingsfilename), ...
            path_regdata, the_inputs.expID_base, the_inputs.expIDs_2align{:});
        
    otherwise
        error('Command not found')
end






end


function [varargout] = check_whitespaces_in_path(varargin)

% Find and correct any path with whitespaces
varargout = varargin;

for iPath=1:nargin
    
    if any(isspace(varargin{iPath}))
        varargout{iPath} = replace(varargin{iPath},' ','\ ');
    end
    
    
end

end