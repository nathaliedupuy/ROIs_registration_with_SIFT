function [REG_RES_n] = edit_rois_regtype(action, REG_RES_n, previous_type, new_type, rois_pair)

switch action
    case 'swap_pair'
        [REG_RES_n] = swap_pair_type(REG_RES_n, previous_type, new_type, rois_pair);
        
    case 'delete_pair'
        [REG_RES_n] = delete_pair_of_rois(REG_RES_n, previous_type, new_type, rois_pair);
        
    case 'add_pair'
        [REG_RES_n] = add_pair_of_rois(REG_RES_n, previous_type, new_type, rois_pair);
        
    otherwise
        error('swap not recognized')
end

end



function [REG_RES_n] = swap_pair_type(REG_RES_n, previous_type, new_type, rois_pair)

% Find the pair
is_rois_pair_2change = REG_RES_n.(previous_type)(1, :) == rois_pair(1);

if any(is_rois_pair_2change)
    % Otherwise it means that it was already moved
    
    % Remove it
    REG_RES_n.(previous_type) = REG_RES_n.(previous_type)(:, ~is_rois_pair_2change);
    
    % Now add it to the new pair type
    REG_RES_n.(new_type) = [REG_RES_n.(new_type), rois_pair];
end
    
end

function [REG_RES_n] = delete_pair_of_rois(REG_RES_n, previous_type, new_type, rois_pair)

if ~isempty(previous_type)
    
    if isempty(REG_RES_n.(previous_type))
        move_pair = false;
    else
        % Find the pair
        is_rois_pair_2change = REG_RES_n.(previous_type)(1, :) == rois_pair(1);
        
        move_pair = any(is_rois_pair_2change);
    end
    
    if move_pair % Otherwise it means that it was already moved
        % Remove it
        REG_RES_n.(previous_type) = REG_RES_n.(previous_type)(:, ~is_rois_pair_2change);
    end
end

if ~isempty(new_type)
    if isempty(REG_RES_n.(new_type){1})
        move_pair = true;
    else
        move_pair = ~any(REG_RES_n.(new_type){1} == rois_pair(1));
    end
    
    
    if move_pair
        REG_RES_n.(new_type){1} = [REG_RES_n.(new_type){1}, rois_pair(1)];
        REG_RES_n.(new_type){2} = [REG_RES_n.(new_type){2}, rois_pair(2)];
    end
end

end

function [REG_RES_n] = add_pair_of_rois(REG_RES_n, previous_type, new_type, rois_pair)

if ~isempty(previous_type)

    if isempty(REG_RES_n.(previous_type){1})
        move_pair = false;
    else
        is_roi_2change = REG_RES_n.(previous_type){1} == rois_pair(1);
        
        move_pair = any(is_roi_2change);
    end
    
    
    if move_pair % Otherwise it means that it was already moved
        REG_RES_n.(previous_type){1} = REG_RES_n.(previous_type){1}(~is_roi_2change);
        
        is_roi_2change = REG_RES_n.(previous_type){2} == rois_pair(2);
        REG_RES_n.(previous_type){2} = REG_RES_n.(previous_type){2}(~is_roi_2change);
    end
end

if ~isempty(new_type)
    if isempty(REG_RES_n.(new_type))
        move_pair = true;
    else
        move_pair = ~any(REG_RES_n.(new_type)(1,:) == rois_pair(1));
    end

    if move_pair
        % Now dd it to the new pair type
        REG_RES_n.(new_type) = [REG_RES_n.(new_type), rois_pair];
    end
end

end