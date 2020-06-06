function [REG_RES, new_flag] = follow_pair_chain(list_sessions, REG_RES, current_flag, rois_pair, idx_3sessions, list_visited_triplet)
% Recursive function
%
% rois_pair: [2 x 1]

all_paired_sessions = {REG_RES.sessions};
% Find the REG_RES corresponding to the new session
idx_regres_1 = find(cellfun(@(x)all(ismember(list_sessions([idx_3sessions(1),idx_3sessions(3)]),x)),all_paired_sessions,'un', true));
idx_regres_2 = find(cellfun(@(x)all(ismember(list_sessions([idx_3sessions(2),idx_3sessions(3)]),x)),all_paired_sessions,'un', true));

[new_flag, roi_match] = find_pairs_in_new_session(list_sessions, REG_RES([idx_regres_1, idx_regres_2]), current_flag, rois_pair, idx_3sessions);

% Add the session triplet to the list
list_visited_triplet = [list_visited_triplet; sort(idx_3sessions, 'ascend')];

new_session_idx = [];
next_idx_3sessions = [];
next_rois_pair = [];

% Look for a new session

% Try to link with first roi from pair
if isempty(roi_match{1})
    if ~isempty(roi_match{2})
        % If there were no match for roi 1, try to link with second roi match
        [new_session_idx] = find_new_connected_session(list_sessions, REG_RES, rois_pair(1), roi_match{2}, list_visited_triplet, idx_3sessions(1), idx_3sessions(3));
        if ~isempty(new_session_idx)
            next_idx_3sessions = [idx_3sessions([1,3]), new_session_idx];
            next_rois_pair = [rois_pair(1); roi_match{2}];
        end
    end
else
    % Try to link with first roi match
    [new_session_idx] = find_new_connected_session(list_sessions, REG_RES, rois_pair(1), roi_match{1}, list_visited_triplet, idx_3sessions(1), idx_3sessions(3));
    if ~isempty(new_session_idx)
        next_idx_3sessions = [idx_3sessions([1,3]), new_session_idx];
        next_rois_pair = [rois_pair(1); roi_match{1}];
    end
end
% If no link found, try to link with second roi from pair
if isempty(new_session_idx) 
    if isempty(roi_match{2})
        if ~isempty(roi_match{1})
            [new_session_idx] = find_new_connected_session(list_sessions, REG_RES, rois_pair(2), roi_match{1}, list_visited_triplet, idx_3sessions(2), idx_3sessions(3));
            if ~isempty(new_session_idx)
                next_idx_3sessions = [idx_3sessions([2,3]), new_session_idx];
                next_rois_pair = [rois_pair(2); roi_match{1}];
            end
        end
    else
        [new_session_idx] = find_new_connected_session(list_sessions, REG_RES, rois_pair(2), roi_match{2}, list_visited_triplet, idx_3sessions(2), idx_3sessions(3));
        if ~isempty(new_session_idx)
            next_idx_3sessions = [idx_3sessions([2,3]), new_session_idx];
            next_rois_pair = [rois_pair(2); roi_match{2}];
        end
    end
end


% If found a new linked session, follow the chain !
if ~isempty(next_idx_3sessions)
    [REG_RES, new_flag] = follow_pair_chain(list_sessions, REG_RES, new_flag, next_rois_pair, next_idx_3sessions, list_visited_triplet);
end

% Once this is done, we can update the REG_RES
new_session = list_sessions{idx_3sessions(3)};
REG_RES(idx_regres_1) = update_regtype_given_newflag(REG_RES(idx_regres_1), rois_pair(1), roi_match{1}, roi_match{2}, new_session, new_flag);
REG_RES(idx_regres_2) = update_regtype_given_newflag(REG_RES(idx_regres_2), rois_pair(2), roi_match{2}, roi_match{1}, new_session, new_flag);


end



function [new_flag, roi_match] = find_pairs_in_new_session(list_sessions, REG_RES_nm, current_flag, rois_pair, idx_3sessions)
% REG_RES_nm : [2 x 1]
% rois_pair : [2 x 1]

new_session = list_sessions{idx_3sessions(3)};
prev_sessions = list_sessions(idx_3sessions([1,2]));

% Find match in the new session for the two rois of the given pair
roi_match = arrayfun(@(rs,roi) find_roi_match(rs.sessions,rs.matched_pairs,new_session,roi), REG_RES_nm, rois_pair, 'un', false);

if current_flag==-1
    % Conflict
    new_flag = current_flag;
else
    % Check new matches
    has_match = [~isempty(roi_match{1}), ~isempty(roi_match{2})];
    if all(~has_match) || (all(has_match) && (roi_match{1}==roi_match{2}))
        % Matched rois are the same
        new_flag = current_flag; % It could be part of complete or incomplete chain
    else
        if all(has_match)
            % The pairs don't match - Conflict
            new_flag = -1;
        else
            % If roiY in session Y had no match in new session, and if roiX
            % from session X was matched to roiN in new session, find out
            % if roiN has a match in session Y
            roi_match_back = find_roi_match(REG_RES_nm(~has_match).sessions, REG_RES_nm(~has_match).matched_pairs, prev_sessions{~has_match}, roi_match{has_match});
            if isempty(roi_match_back)
                % The pair might be missing
                new_flag = 0;
            else
                % The pairs don't match - Conflict
                new_flag = -1;
            end
        end
        
    end
end

end


