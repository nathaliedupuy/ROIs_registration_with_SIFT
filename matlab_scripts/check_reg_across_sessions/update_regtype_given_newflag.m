function [REG_RES_n] = update_regtype_given_newflag(REG_RES_n, roi, roi_match_1, roi_match_2, new_session, new_flag)
% Update the registration of a pair or potential new pair after obtained
% the flag of registration following a chain of pairs.

% We need to prep the pair according to the order in registration.

% Find where the new session is
is_new_session = strcmp(REG_RES_n.sessions, new_session);
% Prep the pair
rois_pair = zeros(2,1);
rois_pair(~is_new_session) = roi;


if new_flag==0 % The chain is missing some pairs, but no conflict
    
    if ~isempty(roi_match_1)
        % Incomplete chain
        rois_pair(is_new_session) = roi_match_1;
        REG_RES_n = edit_rois_regtype('swap_pair', REG_RES_n, 'matched_pairs', 'incomplete_chain', rois_pair);

    elseif ~isempty(roi_match_2)
        % Missed pair
        rois_pair(is_new_session) = roi_match_2;
        
        REG_RES_n = edit_rois_regtype('add_pair', REG_RES_n, 'nonmatched_ROIs', 'missed_pairs', rois_pair);
        
    end
    
else
    % The chain of pairs is either good or conflicting
    if ~isempty(roi_match_1)
        
        rois_pair(is_new_session) = roi_match_1;
        
        if new_flag==1
            REG_RES_n = edit_rois_regtype('swap_pair', REG_RES_n, 'matched_pairs', 'matched_pairs_checked', rois_pair);
        else
            REG_RES_n = edit_rois_regtype('swap_pair', REG_RES_n, 'matched_pairs', 'conflicting_chain', rois_pair);
            
%             % For conflicting we add it to the conflict chain list but do
%             % not remove it from the matched list yet
%             REG_RES_n = edit_rois_regtype('add_pair', REG_RES_n, [], 'conflicting_chain', rois_pair);
        end
        
    end
    
end




end