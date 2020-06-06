function [REG_RES] = check_reg_across_sessions(list_sessions, REG_RES)
% Check if matched pairs are transitive across sessions. Find conflicting
% chains of pairs, and incomplete chains with and missing pairs.

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



% Initialise REG_RES with new fields:
% - a dummy field that will store all valid pairs of matched rois
% - conflicting_chain
% - incomplete_chain
% - missed_pairs
REG_RES(1).matched_pairs_checked = [];
REG_RES(1).conflicting_chain = [];
REG_RES(1).incomplete_chain = [];
REG_RES(1).missed_pairs = [];

for ireg=1:length(REG_RES)
    
    % Get the current sessions indices
    [~, idx_session_pair] = ismember(REG_RES(ireg).sessions, list_sessions);
    
    % Loop over matched rois
    while ~isempty(REG_RES(ireg).matched_pairs)
        
        % Initialise
        current_flag = 1;
        list_visited_triplet = [];
        rois_pair = REG_RES(ireg).matched_pairs(:,1); % [2 x 1]
        
        % Look for a connected session
        [new_session_idx] = find_new_connected_session(list_sessions, REG_RES, rois_pair(1), rois_pair(2), list_visited_triplet, idx_session_pair(1), idx_session_pair(2));
        
        if isempty(new_session_idx)
            % Directly save the pair by removing it from the list of
            % matched pairs and adding it to the temporary field of checked
            % matched pairs.
            REG_RES(ireg) = edit_rois_regtype('swap_pair', REG_RES(ireg), 'matched_pairs', 'matched_pairs_checked', rois_pair);

        else
            [REG_RES, new_flag] = follow_pair_chain(list_sessions, REG_RES, current_flag, rois_pair, [idx_session_pair, new_session_idx], list_visited_triplet);            
            
            % Final update with the initial pair
            REG_RES(ireg) = update_regtype_given_newflag(REG_RES(ireg), rois_pair(1), rois_pair(2), [], REG_RES(ireg).sessions{2}, new_flag);


        end
        
    end
    
    
end


for ireg=1:length(REG_RES)
    REG_RES(ireg).matched_pairs = REG_RES(ireg).matched_pairs_checked; % Only keep the ones that checked out
end
% Finally, remove the temporary field
REG_RES = rmfield(REG_RES,'matched_pairs_checked');

end




