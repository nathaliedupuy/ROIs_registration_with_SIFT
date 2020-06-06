function [idx_new_session, roi_match, new_reg_idx] = find_new_connected_session(list_sessions, REG_RES, roi_1, roi_2, list_visited_triplet, idx_session_1, idx_session_2)
% Find a new connected session to a given session or a pair of sessions. A
% connected session is a session with a roi match.
%
% - For a given roi in a given session (leave empty the other roi and
%   session inputs), find a new session with a roi match.
%
% - For a given pair of rois from two session, find a new session with a
%   roi match to one of the input rois. the roi match has to be valid in
%   the registration with both original sessions (i.e. not deleted not out
%   of fov).

idx_session_pair = [idx_session_1, idx_session_2];

num_sessions = length(list_sessions);
list_ses_idx = setdiff(1:num_sessions, idx_session_pair);
num_visited = size(list_visited_triplet,1);

all_paired_session = {REG_RES.sessions};

idx_new_session = [];
roi_match = [];
new_reg_idx = [];
for idxSes=list_ses_idx
    
    triplet_was_explored = ~isempty(list_visited_triplet) && any(all(repmat(sort([idxSes,idx_session_pair],'ascend'), [num_visited,1]) == list_visited_triplet ,2));
    if ~triplet_was_explored
        % Check if any of the rois is connected to this new session
        is_valid = true;
        
        % Try with first roi
        if ~isempty(roi_1)
            is_regres1 = cellfun(@(x)all(ismember(list_sessions([idx_session_pair(1),idxSes]),x)),all_paired_session,'un', true);
            is_regres = is_regres1;
            [roi_match] = find_roi_match(REG_RES(is_regres1).sessions, REG_RES(is_regres1).matched_pairs, list_sessions(idxSes), roi_1);
        end
        
        if ~isempty(roi_2)
            is_regres2 = cellfun(@(x)all(ismember(list_sessions([idx_session_pair(2),idxSes]),x)),all_paired_session,'un', true);
            if ~isempty(roi_match)
                % We have to make sure that the roi match was not deleted or out of
                % fov in the other registration
                is_valid = check_roi_is_valid(REG_RES(is_regres2).sessions, list_sessions(idxSes), roi_match, ...
                    REG_RES(is_regres2).outofFOV_ROIs, REG_RES(is_regres2).deleted_ROIs);
            else
                % Try with other session
                is_regres = is_regres2;
                [roi_match] = find_roi_match(REG_RES(is_regres2).sessions, REG_RES(is_regres2).matched_pairs, list_sessions(idxSes), roi_2);
                if ~isempty(roi_match) && ~isempty(roi_1)
                    % We have to make sure that the roi match was not deleted or out of
                    % fov in the other registration
                    is_valid = check_roi_is_valid(REG_RES(is_regres1).sessions, list_sessions(idxSes), roi_match, ...
                        REG_RES(is_regres1).outofFOV_ROIs, REG_RES(is_regres1).deleted_ROIs);
                end
            end
        end
       
        if ~isempty(roi_match) && is_valid
            idx_new_session = idxSes;
            new_reg_idx = find(is_regres);
            break
        end
    end
end

end

function [is_valid] = check_roi_is_valid(sessions, new_session, roi, outofFOV_ROIs, deleted_ROIs)

is_new_session = strcmp(sessions, new_session);

is_valid = true;
if ~isempty(outofFOV_ROIs{is_new_session})
    is_valid = ~any(roi==outofFOV_ROIs{is_new_session});
end
if is_valid && ~isempty(deleted_ROIs{is_new_session})
    is_valid = ~any(roi==deleted_ROIs{is_new_session});
end

end
