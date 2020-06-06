function [roi_match] = find_roi_match(sessions, matched_pairs, new_session, roi)

is_new_session = strcmp(sessions, new_session);
idx_match = find(roi==matched_pairs(~is_new_session,:));
if isempty(idx_match)
    roi_match = [];
else
    roi_match = matched_pairs(is_new_session,idx_match);
end

end