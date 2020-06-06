function [REG_RES] = gui_update_REG_RES_across_session(list_sessions, REG_RES, current_registration, idx_session_pair, rois, inspect_type, varargin)
% Check if match are trasitive across sessions, and flag conflicts.

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

num_ses = length(list_sessions);

% Make sure rois is [2x1]
rois = rois(:);

switch inspect_type
    
    case 'matched_pairs' % -- Delete a matched pair
        
        % Remove pair from 'matched_pairs' but do not add the rois to
        % nonmatched yet
        REG_RES(current_registration) = edit_rois_regtype('delete_pair', REG_RES(current_registration), 'matched_pairs', [], rois);

        if num_ses==2
            % If 2 sessions, directly move the rois to 'nonmatched_ROIs'
            REG_RES(current_registration) = edit_rois_regtype('delete_pair', REG_RES(current_registration), [], 'nonmatched_ROIs', rois);
            
        else
            % If more than 2 sessions, we might have to do a recursive search
            % to edit all connected pairs
            
            % The search will be only through the matched_pairs
            [dummyregres] = initialise_dummyregres({REG_RES.sessions}, {REG_RES.matched_pairs}, {REG_RES.nonmatched_ROIs}, {REG_RES.outofFOV_ROIs}, {REG_RES.deleted_ROIs});
            
            % Look for a connected session
            [new_session_idx] = find_new_connected_session(list_sessions, dummyregres, rois(1), rois(2), [], idx_session_pair(1), idx_session_pair(2));
            
            if isempty(new_session_idx)
                % If did not find a connection to other sessions, directly move
                % the rois to 'nonmatched_ROIs'
                REG_RES(current_registration) = edit_rois_regtype('delete_pair', REG_RES(current_registration), [], 'nonmatched_ROIs', rois);
                
            else
                % If there is a connection, we need to ID all connected
                % pairs and move them to incomplete chain.
                current_flag = 1;
                [dummyregres, newflag] = follow_pair_chain(list_sessions, dummyregres, current_flag, rois, [idx_session_pair, new_session_idx], []);
                % Note: newflag should be = 1
                if newflag==-1
                    error('something went wrong')
                end
                for ireg=1:length(REG_RES)
                    if ireg==current_registration
                        % We add the newly deleted pair to the missings...
                        REG_RES(current_registration).missed_pairs = [REG_RES(current_registration).missed_pairs, rois]; 
                    else
                        % We move the matched_pairs_checked which are
                        % related to the newly deleted pair to the
                        % incomplete_chain
                        REG_RES(ireg).incomplete_chain = [REG_RES(ireg).incomplete_chain, dummyregres(ireg).matched_pairs_checked];
                        % Update the matched pairs - Should only remain the
                        % ones that were not visited during the search
                        REG_RES(ireg).matched_pairs = dummyregres(ireg).matched_pairs; 
                    end
                end
            end
        end
        
        
        
    case 'conflicting_chain' % -- Delete a matched pair part of a conflicting chain of pairs
        
        % Remove pair from 'conflicting_chain' and Directly move the rois
        % to 'nonmatched_ROIs'
        REG_RES(current_registration) = edit_rois_regtype('delete_pair', REG_RES(current_registration), 'conflicting_chain', 'nonmatched_ROIs', rois);

        % Note, by default there are more than 2 sessions to have conflict.
        % Also, there should be at least one connected session; otherwise
        % this match wouldn't be in conflicting chain of pairs.
        
        % Start a dummy registration to see if removing this pair resolved
        % anything; the search will be only through the conflicting_chain.
        
        % We first check the chain starting from first roi of the deleted
        % pair, and update REG_RES.
        [REG_RES] = search_update_conflict_chain(list_sessions, REG_RES, rois(1), idx_session_pair(1));
        
        % We then check starting from second roi of the deleted pair.
        [REG_RES] = search_update_conflict_chain(list_sessions, REG_RES, rois(2), idx_session_pair(2));
        
        
    case 'incomplete_chain' % -- Delete a matched pair part of a incomplete chain of pairs
        
        % Note, by default there are more than 2 sessions to have missing
        % pairs. Also, there should be at least one connected session;
        % otherwise this match wouldn't be in an incomplete chain of pairs.
        
        % First run a recursive search to identify the pool of missing
        % pairs that are related to it.
        
        POOL_nonmatched_ROIs = cell(length(REG_RES),1); % has same format as nonmatched_ROIs
        for ireg=1:length(REG_RES)
            POOL_nonmatched_ROIs{ireg} = cell(1,2);
            if ~isempty(REG_RES(ireg).missed_pairs)
                POOL_nonmatched_ROIs{ireg}{1} = REG_RES(ireg).missed_pairs(1,:);
                POOL_nonmatched_ROIs{ireg}{2} = REG_RES(ireg).missed_pairs(2,:);
            end 
        end
        
        [dummyregres] = initialise_dummyregres({REG_RES.sessions},{REG_RES.incomplete_chain}, POOL_nonmatched_ROIs, {REG_RES.outofFOV_ROIs}, {REG_RES.deleted_ROIs});
        % Look for a connected session
        [new_session_idx] = find_new_connected_session(list_sessions, dummyregres, rois(1), rois(2), [], idx_session_pair(1), idx_session_pair(2));
        % Run the search.
        current_flag = 1;
        [dummyregres, ~] = follow_pair_chain(list_sessions, dummyregres, current_flag, rois, [idx_session_pair, new_session_idx], []);
        
        % Now use the results to identify the pool of missed_pairs related
        % to the deleted pair, which we will use as nonmatched.
        POOL_nonmatched_ROIs = cell(length(REG_RES),1); % has same format as nonmatched_ROIs
        for ireg=1:length(REG_RES)
            if ireg==current_registration
                % We add the deleted pair to the pool
                POOL_nonmatched_ROIs{ireg} = {rois(1), rois(2)};
            else
                if isempty(REG_RES(ireg).missed_pairs) || isempty(dummyregres(ireg).missed_pairs)
                    POOL_nonmatched_ROIs{ireg} = {[], []};
                else
                    % We remove from REG_RES the missed pairs that were
                    % identified as part of the chain.
                    is_in_pool = REG_RES(ireg).missed_pairs(1,:) == dummyregres(ireg).missed_pairs(1);
                    REG_RES(ireg).missed_pairs = REG_RES(ireg).missed_pairs(:,~is_in_pool);
                    POOL_nonmatched_ROIs{ireg} = {dummyregres(ireg).missed_pairs(1), dummyregres(ireg).missed_pairs(2)};
                end
            end
        end
        
        % Remove pair from 'incomplete_chain' but do not add the rois to
        % nonmatched yet
        REG_RES(current_registration) = edit_rois_regtype('delete_pair', REG_RES(current_registration), 'incomplete_chain', [], rois);
                
        % Start a dummy registration to see if removing this pair changed
        % anything; the search will be only through the incomplete_chain.
        % We use the POOL_missed_pairs as nonmatched_ROIs for the search.
        
        % We first check the chain starting from the first roi of the
        % deleted pair, update REG_RES, and get the new nonmatched pool for
        % the next round.
        [REG_RES, POOL_nonmatched_ROIs] = search_update_incomplete_chain(list_sessions, REG_RES, POOL_nonmatched_ROIs, rois(1), idx_session_pair(1));
        
        % Next we repeat with roi of other session
        [REG_RES, POOL_nonmatched_ROIs] = search_update_incomplete_chain(list_sessions, REG_RES, POOL_nonmatched_ROIs, rois(2), idx_session_pair(2));
        
        % Finally, we need to update the nonmatched_ROIs with the
        % missed_pairs that didn't make the cut anywhere.
        for ireg=1:length(REG_RES)
            REG_RES(ireg).nonmatched_ROIs{1} = [REG_RES(ireg).nonmatched_ROIs{1}, POOL_nonmatched_ROIs{ireg}{1}];
            REG_RES(ireg).nonmatched_ROIs{2} = [REG_RES(ireg).nonmatched_ROIs{2}, POOL_nonmatched_ROIs{ireg}{2}];
        end
        
        
    case 'missed_pairs' 
        % --  Couple of options here
        % --  (1) Accept a matched missing pair part of a incomplete chain of pairs
        % --  (2) Delete roi
        % --  (3) Out of fov roi
        
        % Note, by default there are more than 2 sessions to have missing
        % pairs. Also, there should be at least one connected session;
        % otherwise this match wouldn't be in an incomplete chain of pairs.
        
        if isempty(varargin)
            % Remove pair from 'missed_pairs' and add to incomplete_chain
            REG_RES(current_registration) = edit_rois_regtype('swap_pair', REG_RES(current_registration), 'missed_pairs', 'incomplete_chain', rois);
        else
            % Remove pair from 'missed_pairs' and delete both rois
            REG_RES(current_registration) = edit_rois_regtype('delete_pair', REG_RES(current_registration), 'missed_pairs', varargin{1}, rois);
        end
        
        % Merge nonmatched_ROIs and missed pairs to check if completed
        % chain or discovered new missing pairs
        POOL_nonmatched_ROIs = cell(length(REG_RES),1); % has same format as nonmatched_ROIs
        % And remove them.
        for ireg=1:length(REG_RES)
            POOL_nonmatched_ROIs{ireg} = REG_RES(ireg).nonmatched_ROIs;
            if ~isempty(REG_RES(ireg).missed_pairs)
                POOL_nonmatched_ROIs{ireg}{1} = [POOL_nonmatched_ROIs{ireg}{1}, REG_RES(ireg).missed_pairs(1,:)];
                POOL_nonmatched_ROIs{ireg}{2} = [POOL_nonmatched_ROIs{ireg}{2}, REG_RES(ireg).missed_pairs(2,:)];
            end 
        end
        
        % Run a recursive to check if the chain is now complete
        [dummyregres] = initialise_dummyregres({REG_RES.sessions},{REG_RES.incomplete_chain}, POOL_nonmatched_ROIs, {REG_RES.outofFOV_ROIs}, {REG_RES.deleted_ROIs});

        % Look for a connected session
        [new_session_idx] = find_new_connected_session(list_sessions, dummyregres, rois(1), rois(2), [], idx_session_pair(1), idx_session_pair(2));
        % Run the search.
        current_flag = 1;
        [dummyregres, new_flag] = follow_pair_chain(list_sessions, dummyregres, current_flag, rois, [idx_session_pair, new_session_idx], []);
        
        % Now we need to change the pairs of incomplete chains that have been
        % moved to matched_pairs. We need to do this only if the result
        % 'new_flag' has changed from 0->1
        if new_flag==1
            for ireg=1:length(REG_RES)
                if ireg==current_registration
                    % If it's a match, we move it to the right place
                    REG_RES(current_registration) = edit_rois_regtype('swap_pair', REG_RES(current_registration), 'incomplete_chain', 'matched_pairs', rois);
                else
                    if ~isempty(dummyregres(ireg).matched_pairs_checked)
                        REG_RES(ireg) = edit_rois_regtype('swap_pair', REG_RES(ireg), 'incomplete_chain', 'matched_pairs', dummyregres(ireg).matched_pairs_checked);
                    end
                end
            end
        else
            % We might have discovered new missed pairs
            for ireg=1:length(REG_RES)
                if ireg==current_registration
                    continue
                else
                    if ~isempty(dummyregres(ireg).missed_pairs)
                        % Remove from nonmatched_ROIs and add to
                        % missed_pairs
                        REG_RES(ireg) = edit_rois_regtype('add_pair', REG_RES(ireg), 'nonmatched_ROIs', 'missed_pairs', dummyregres(ireg).missed_pairs);                        
                    end
                end
            end
        end
        
        
    case 'nonmatched_ROIs' % -- Match 2 nonmatched ROIs
        
        % Remove rois from nonmatched
        REG_RES(current_registration) = edit_rois_regtype('add_pair', REG_RES(current_registration), 'nonmatched_ROIs', [], rois);  

        if num_ses==2
            % If only 2 sessions, directly move the rois to
            % 'matched_pairs'
            REG_RES(current_registration).matched_pairs = [REG_RES(current_registration).matched_pairs, rois]; % Add the new pair
            
        else
            
            % If more than 2 sessions, we might have to do a recursive search
            % for connected pairs
            
            % Given the new pair, there are couple of options:
            % - it is a simple match
            % - it is linked to matched_pairs, incomplete_chains, or
            % conflicting_chains; and on top of that it can link different
            % types of pairs... => 2 case scenarios: it leads to
            % incomplete_chain or conflicting_chain
            
            % Run a global search - prep the pools
            POOL_matched_pairs = cell(length(REG_RES),1); 
            POOL_nonmatched_ROIs = cell(length(REG_RES),1); % has same format as nonmatched_ROIs
            for ireg=1:length(REG_RES)
                POOL_nonmatched_ROIs{ireg} = REG_RES(ireg).nonmatched_ROIs;
                if ~isempty(REG_RES(ireg).missed_pairs)
                    POOL_nonmatched_ROIs{ireg}{1} = [POOL_nonmatched_ROIs{ireg}{1}, REG_RES(ireg).missed_pairs(1,:)];
                    POOL_nonmatched_ROIs{ireg}{2} = [POOL_nonmatched_ROIs{ireg}{2}, REG_RES(ireg).missed_pairs(2,:)];
                end
                if ireg==current_registration
                    % We add the new pair to the pool
                    POOL_matched_pairs{ireg} = [REG_RES(ireg).matched_pairs, rois, REG_RES(ireg).conflicting_chain, REG_RES(ireg).incomplete_chain];
                else
                    POOL_matched_pairs{ireg} = [REG_RES(ireg).matched_pairs, REG_RES(ireg).conflicting_chain, REG_RES(ireg).incomplete_chain];
                    
                end
            end
            
            [dummyregres] = initialise_dummyregres({REG_RES.sessions},POOL_matched_pairs, POOL_nonmatched_ROIs, {REG_RES.outofFOV_ROIs}, {REG_RES.deleted_ROIs});
            
            % Look for a connected session
            [new_session_idx] = find_new_connected_session(list_sessions, dummyregres, rois(1), rois(2), [], idx_session_pair(1), idx_session_pair(2));
            
            if isempty(new_session_idx)
                % If did not find a connection to other sessions, directly
                % move the rois to 'matched_pairs'
                REG_RES(current_registration).matched_pairs = [REG_RES(current_registration).matched_pairs, rois]; % Add the new pair
                
            else
                % Otherwise we have to run recursive search to see what has
                % been changed by the new pair.
                current_flag = 1;
                [dummyregres, new_flag] = follow_pair_chain(list_sessions, dummyregres, current_flag, rois, [idx_session_pair, new_session_idx], []);
                
                % Note: Because the pair is linked, new_flag should be
                % either 0 ot -1, and matched_pairs_checked in dummyregres
                % should be empty.
                if new_flag==-1
                    f={'conflicting_chain'}; f0=f{1};
                else
                    f={'incomplete_chain', 'missed_pairs'};
                end
                
                % One last edit for the pair in case was not done during
                % recursive search
                dummyregres(current_registration) = edit_rois_regtype('swap_pair', dummyregres(current_registration), 'matched_pairs', f{1}, dummyregres(current_registration).(f{1}));

                % Now loop to update with the changes
                for ireg=1:length(REG_RES)
                    % We look for new or moved pairs
                    if new_flag==0
                        b = [~isempty(dummyregres(ireg).(f{1})), ~isempty(dummyregres(ireg).(f{2}))];
                        if any(b)
                            f0 = f{b};
                        else
                            f0 = f{1};
                        end
                    end
                    if ~isempty(dummyregres(ireg).(f0))
                        % We have to search for the original location
                        % of the pair
                        [type_origin] = look_for_rois_pair_origin(REG_RES(ireg), dummyregres(ireg).(f0));
                        if strcmp(type_origin, 'nonmatched_ROIs')
                            % Remove them from nonmatched and Add to new pair type
                            REG_RES(ireg) = edit_rois_regtype('add_pair', REG_RES(ireg), 'nonmatched_ROIs', f0, dummyregres(ireg).(f0));
                        else
                            % Swap, but only if was changed...
                            if ~strcmp(type_origin, f0)
                                REG_RES(ireg) = edit_rois_regtype('swap_pair', REG_RES(ireg), type_origin, f0, dummyregres(ireg).(f0));
                            end
                        end
                    end
                end
            end

        end
        
        
    otherwise % -- Add a new roi, either from outofFOV_ROIs or deleted_ROIs
        
        % Get session and roi
        is_ses = [~isempty(rois{1}), ~isempty(rois{2})];
        sesname = REG_RES(current_registration).sessions{is_ses};
        
        % Move the new roi to nonmatched_ROIs
        [REG_RES(current_registration)] = add_roi(REG_RES(current_registration), rois{is_ses}, is_ses, inspect_type);
                
        if num_ses>2
            % Prep the pools for global search
            POOL_matched_pairs = cell(length(REG_RES),1);
            POOL_nonmatched_ROIs = cell(length(REG_RES),1); % has same format as nonmatched_ROIs
            for ireg=1:length(REG_RES)
                POOL_matched_pairs{ireg} = [REG_RES(ireg).matched_pairs, REG_RES(ireg).conflicting_chain, REG_RES(ireg).incomplete_chain];
                POOL_nonmatched_ROIs{ireg} = REG_RES(ireg).nonmatched_ROIs;
                if ~isempty(REG_RES(ireg).missed_pairs)
                    POOL_nonmatched_ROIs{ireg}{1} = [POOL_nonmatched_ROIs{ireg}{1}, REG_RES(ireg).missed_pairs(1,:)];
                    POOL_nonmatched_ROIs{ireg}{2} = [POOL_nonmatched_ROIs{ireg}{2}, REG_RES(ireg).missed_pairs(2,:)];
                end
            end
            [dummyregres] = initialise_dummyregres({REG_RES.sessions},POOL_matched_pairs, POOL_nonmatched_ROIs, {REG_RES.outofFOV_ROIs}, {REG_RES.deleted_ROIs});
            
            % We seach for any connection to this roi, and run a search if
            % found one.
            [dummyregres, check_if_changed] = search_update_add_roi(list_sessions, dummyregres, current_registration, sesname, rois{is_ses});

            if check_if_changed
                f={'matched_pairs_checked','conflicting_chain','incomplete_chain', 'missed_pairs'};
                % Now loop to update with the changes
                for ireg=1:length(REG_RES)
                    % We look changes - note that there should only be one
                    % change per registration, and only affect pairs.
                    b = cellfun(@(s)~isempty(dummyregres(ireg).(s)), f , 'un', true);
                    f0 = f{b}; 
                    if ~isempty(f0)
                        f00=f0;
                        if strcmp(f00, 'matched_pairs_checked')
                            f00 = 'matched_pairs';
                        end
                        % We have to search for the original location
                        % of the pair
                        [type_origin] = look_for_rois_pair_origin(REG_RES(ireg), dummyregres(ireg).(f0));
                        if strcmp(type_origin, 'nonmatched_ROIs')
                            % Remove them from nonmatched and Add to new pair type
                            REG_RES(ireg) = edit_rois_regtype('add_pair', REG_RES(ireg), 'nonmatched_ROIs', f00, dummyregres(ireg).(f0));                           
                        else
                            % Swap, but only if was changed...
                            if ~strcmp(type_origin, f00)
                                REG_RES(ireg) = edit_rois_regtype('swap_pair', REG_RES(ireg), type_origin, f00, dummyregres(ireg).(f0));
                            end
                        end                        
                    end
                end
            end
        end

end

end



function [REG_RES, check_if_changed] = search_update_add_roi(list_sessions, REG_RES, current_registration, sesname, roi)

check_if_changed = false;

for ireg=1:length(REG_RES)
    
    if ireg==current_registration
        % We know there are no pair here, since we just added the new roi.
        continue
    end
    % Check has the session 
    [b, k_ses] = ismember(sesname, REG_RES(ireg).sessions);
    if ~b
        continue
    end

    % See if find any with the new roi
    has_pair_with_the_roi = REG_RES(ireg).matched_pairs(k_ses,:) == roi;
    if ~any(has_pair_with_the_roi)
        continue
    end
    
    
    % Take the pair:
    rois_pair = REG_RES(ireg).matched_pairs(:, has_pair_with_the_roi);
    
    % Get the current sessions indices
    [~, idx_session_pair] = ismember(REG_RES(ireg).sessions, list_sessions);
    
    % Run search
    current_flag = 1;
    list_visited_triplet = [];
    
    % Look for a connected session
    [new_session_idx] = find_new_connected_session(list_sessions, REG_RES, rois_pair(1), rois_pair(2), list_visited_triplet, idx_session_pair(1), idx_session_pair(2));
    
    if isempty(new_session_idx)
        % Directly save the pair by removing it from the list of matched
        % pairs and adding it to the temporary field of checked matched
        % pairs.
        REG_RES(ireg) = edit_rois_regtype('swap_pair', REG_RES(ireg), 'matched_pairs', 'matched_pairs_checked', rois_pair);
        
    else
        if ~check_if_changed
            check_if_changed = true;
        end
        
        [REG_RES, new_flag] = follow_pair_chain(list_sessions, REG_RES, current_flag, rois_pair, [idx_session_pair, new_session_idx], list_visited_triplet);
        
        % Final update with the initial pair
        REG_RES(ireg) = update_regtype_given_newflag(REG_RES(ireg), rois_pair(1), rois_pair(2), [], REG_RES(ireg).sessions{2}, new_flag);
        
    end

end

end

function [REG_RES, POOL_nonmatched_ROIs] = search_update_incomplete_chain(list_sessions, REG_RES, POOL_nonmatched_ROIs, roi_root, idx_session_root)

% Start a dummy registration to see if removing this pair changed anything;
% the search will be only through the incomplete_chain. We use the
% POOL_missed_pairs as nonmatched_ROIs for the search.
[dummyregres] = initialise_dummyregres({REG_RES.sessions},{REG_RES.incomplete_chain}, POOL_nonmatched_ROIs, {REG_RES.outofFOV_ROIs}, {REG_RES.deleted_ROIs});

% Look for a connected session from one of the rois.
[test_session_idx, roi_match, test_reg_idx] = find_new_connected_session(list_sessions, dummyregres, roi_root, [], [], idx_session_root, []);

if ~isempty(test_session_idx)
    % Look again for a connected session
    [new_session_idx] = find_new_connected_session(list_sessions, dummyregres, roi_root, roi_match, [], idx_session_root, test_session_idx);
    
    if isempty(new_session_idx)
        % We can accept the pair as a match if no connection was found
        rois_pair = get_roi_pair(REG_RES(test_reg_idx), 'incomplete_chain', list_sessions{idx_session_root}, roi_root);
        REG_RES(test_reg_idx) = edit_rois_regtype('swap_pair', REG_RES(test_reg_idx), 'incomplete_chain', 'matched_pairs', rois_pair);
        
    else
        % Otherwise we have to search recursively
        current_flag = 1;
        rois_pair = [roi_root; roi_match];
        [dummyregres, new_flag] = follow_pair_chain(list_sessions, dummyregres, current_flag, rois_pair, [idx_session_root, test_session_idx, new_session_idx], []);
        
        % Should not have generated conflict
        if new_flag==-1 
            error('Problem registration update incomplete chain')
        end
        % Last edit, in case was not done during recursive search
        % Final update with the initial pair
        dummyregres(test_reg_idx) = update_regtype_given_newflag(dummyregres(test_reg_idx), rois_pair(1), rois_pair(2), [], list_sessions{test_session_idx}, new_flag);
        
        % Now we need to merge the results
        for ireg=1:length(REG_RES)
            % We might have new matched_pairs
            REG_RES(ireg).matched_pairs = [REG_RES(ireg).matched_pairs, dummyregres(ireg).matched_pairs_checked];
            % For the incomplete chain, we have the ones that are still
            % incomplete after the dummy registration, and we have to merge
            % back the ones that we used to init the 'matched_pairs' in the
            % dummy and were not visited.
            REG_RES(ireg).incomplete_chain = [dummyregres(ireg).matched_pairs, dummyregres(ireg).incomplete_chain];
            % We save the missed_pairs, the missing links for the
            % remaining incomplete chains.
            REG_RES(ireg).missed_pairs = [REG_RES(ireg).missed_pairs, dummyregres(ireg).missed_pairs];
            
            % And finally we save the nonmatched as the new POOL
            POOL_nonmatched_ROIs{ireg} = dummyregres(ireg).nonmatched_ROIs;
            
        end
    end
end


end

function [REG_RES] = search_update_conflict_chain(list_sessions, REG_RES, roi_root, idx_session_root)

% Start a dummy registration to see if removing this pair resolved
% anything; the search will be only through the conflicting_chain.
[dummyregres] = initialise_dummyregres({REG_RES.sessions},{REG_RES.conflicting_chain}, {REG_RES.nonmatched_ROIs}, {REG_RES.outofFOV_ROIs}, {REG_RES.deleted_ROIs});

% Look for a connected session from one of the rois.
[test_session_idx, roi_match, test_reg_idx] = find_new_connected_session(list_sessions, dummyregres, roi_root, [], [], idx_session_root, []);

if ~isempty(test_session_idx)
    
    % Look again for a connected session
    [new_session_idx] = find_new_connected_session(list_sessions, dummyregres, roi_root, roi_match, [], idx_session_root, test_session_idx);
    
    if isempty(new_session_idx)
        % If the pair is not connected to anything else, we can accept the
        % pair as a match
        rois_pair = get_roi_pair(REG_RES(test_reg_idx), 'conflicting_chain', list_sessions{idx_session_root}, roi_root);
        REG_RES(test_reg_idx) = edit_rois_regtype('swap_pair', REG_RES(test_reg_idx), 'conflicting_chain', 'matched_pairs', rois_pair);
        
    else
        % Otherwise we have to search recursively
        current_flag = 1;
        rois_pair = [roi_root; roi_match];
        [dummyregres, new_flag] = follow_pair_chain(list_sessions, dummyregres, current_flag, rois_pair, [idx_session_root, test_session_idx, new_session_idx], []);
        
        % Last edit, in case was not done during recursive search
        % Final update with the initial pair
        dummyregres(test_reg_idx) = update_regtype_given_newflag(dummyregres(test_reg_idx), rois_pair(1), rois_pair(2), [], list_sessions{test_session_idx}, new_flag);
                
        % Now we need to merge the results
        for ireg=1:length(REG_RES)
            % In case we have resolved conflicts, we can add the newly
            % matched pairs / incomplete / missed
            REG_RES(ireg).matched_pairs = [REG_RES(ireg).matched_pairs, dummyregres(ireg).matched_pairs_checked];
            REG_RES(ireg).incomplete_chain = [REG_RES(ireg).incomplete_chain, dummyregres(ireg).incomplete_chain];
            REG_RES(ireg).missed_pairs = [REG_RES(ireg).missed_pairs, dummyregres(ireg).missed_pairs];
            % If new missing pairs have been formed, these are no
            % longer in nonmatched_ROIs, so we need to update.
            REG_RES(ireg).nonmatched_ROIs = dummyregres(ireg).nonmatched_ROIs;
            % For the conflicting, we have the ones that are still
            % conflicting after the dummy registration, and we have to
            % merge back the ones that we used to init the
            % 'matched_pairs' in the dummy and were not visited.
            REG_RES(ireg).conflicting_chain = [dummyregres(ireg).matched_pairs, dummyregres(ireg).conflicting_chain];
        end
    end
    
end


end


function [rois_pair] = get_roi_pair(REG_RES_n, regtype, session_A, roi_A)

% ID order of session
is_session_A = strcmp(REG_RES_n.sessions, session_A);
% Find the pair
is_rois_pair_2change = REG_RES_n.(regtype)(is_session_A, :) == roi_A;
% Save the pair in correct order of registration
rois_pair = REG_RES_n.(regtype)(:, is_rois_pair_2change);

end

function [REG_RES_n] = add_roi(REG_RES_n, roi, is_ses, type_origin)

is_roi_2change = REG_RES_n.(type_origin){is_ses} == roi;

REG_RES_n.(type_origin){is_ses} = REG_RES_n.(type_origin){is_ses}(~is_roi_2change);
REG_RES_n.nonmatched_ROIs{is_ses} = ...
    [REG_RES_n.nonmatched_ROIs{is_ses}, roi];

end

function [type_origin] = look_for_rois_pair_origin(REG_RES_n, tracked_rois)

Slist = {'matched_pairs','conflicting_chain','incomplete_chain','missed_pairs'};

% Find the roi in the original registration
roi__found = false; k=0;
while ~roi__found && k<4
    k=k+1;
    if ~isempty(REG_RES_n.(Slist{k}))
        roi__found = any(tracked_rois(1) == REG_RES_n.(Slist{k})(1,:));
        if roi__found
            % Get the new popup menu for type inspection
            type_origin = Slist{k};
        end
    end
end
% otherwise It has to be in nonmatched
if ~roi__found
    type_origin = 'nonmatched_ROIs';
end

end

function [dummyregres] = initialise_dummyregres(sessions, matched_pairs, nonmatched_ROIs, outofFOV_ROIs, deleted_ROIs)

% Initialise a dummy REG_RES with new fields:
% - a field with some matched_pairs of interest
% - a field with nonmatched_ROIs
% Then fields that will be filled during recursive search
% - a field that will store all valid pairs of matched rois
% - conflicting_chain
% - incomplete_chain
% - missed_pairs

% Note: we might concat several datasets (e.g. matched, conflict, etc...)

n_reg = length(matched_pairs);

sessions = sessions(:);
matched_pairs = matched_pairs(:);
nonmatched_ROIs = nonmatched_ROIs(:);
outofFOV_ROIs = outofFOV_ROIs(:);
deleted_ROIs = deleted_ROIs(:);

dummyregres = struct('sessions', sessions, 'matched_pairs', matched_pairs, 'nonmatched_ROIs', nonmatched_ROIs, ...
    'matched_pairs_checked', cell(n_reg,1), 'conflicting_chain', cell(n_reg,1), ...
    'incomplete_chain', cell(n_reg,1), 'missed_pairs', cell(n_reg,1), ...
    'outofFOV_ROIs', outofFOV_ROIs, 'deleted_ROIs', deleted_ROIs);
end

