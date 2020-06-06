function [BW] = create_roisbinarymask(ROIs_XY, d1, d2, out_as_sparse)
% d1: rows (y)
% d2: columns (x)

if iscell(ROIs_XY)
    BW = cellfun(@(coord_xy)poly2mask(coord_xy(:,1),coord_xy(:,2),d1,d2), ROIs_XY, 'UniformOutput', false);
    
    BW = cat(3, BW{:});
    BW = reshape(BW, [d1*d2, length(ROIs_XY)]);
    
else
    
    BW = poly2mask(ROIs_XY(:,1),ROIs_XY(:,2),d1,d2);
    BW = reshape(BW, [d1*d2, 1]);
end

if nargin ==3 || out_as_sparse
    BW = sparse(BW);
end



end
