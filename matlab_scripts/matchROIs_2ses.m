function [HMatch_results] = matchROIs_2ses(aligned_ROIs_XY, ROIsIDs, imgsize, options)

% NOTE: PART OF CODE ADAPTED FROM CAIMAN PIPELINE 'REGISTER_ROIS.M'


d1 = imgsize(1);
d2 = imgsize(2);

% -------------------------------------------------------------------------
% -- OPTIONS
% -------------------------------------------------------------------------

if ~isfield(options,'dist_exp') || isempty(options.dist_exp); options.dist_exp = 1; end
if ~isfield(options,'dist_thr') || isempty(options.dist_thr); options.dist_thr = 0.5; end
if ~isfield(options,'dist_overlap_thr') || isempty(options.dist_overlap_thr); options.dist_overlap_thr = 0.8; end

% SEE in relevant section for MC options


% -------------------------------------------------------------------------

% Create the sparse, binary masks
M = cell(1,2);
for iSes = 1:2
    
    M{iSes} = create_roisbinarymask(aligned_ROIs_XY{iSes}, d1, d2); % sparse [d x nRois_1]
    
end


num_rois_1 = length(aligned_ROIs_XY{1});
num_rois_2 = length(aligned_ROIs_XY{2});

% Determine distance matrix between M1 and M2
D = zeros(num_rois_1, num_rois_2);
for i = 1:num_rois_1
    for j = 1:num_rois_2
        
        overlap = nnz( M{1}(:,i) & M{2}(:,j) );
        totalarea = nnz( M{1}(:,i) | M{2}(:,j) );
        smallestROI = min(nnz(M{1}(:,i)),nnz(M{2}(:,j)));
        
        D(i,j) = 1 - (overlap/totalarea)^options.dist_exp;
        
        if overlap >= options.dist_overlap_thr*smallestROI
            D(i,j) = 0;
        end
        
    end
end
HMatch_results.score = D; % save it for later!
D(D>options.dist_thr) = Inf;

% Apply Hungarian algo to optimize ROIs matching
R = Hungarian(D);

[match_1,match_2] = find(R);

if isempty(ROIsIDs)
    HMatch_results.matched_pairs = [match_1' ; match_2'];
    HMatch_results.nonmatched{1,1} = setdiff(1:num_rois_1, match_1);
    HMatch_results.nonmatched{1,2} = setdiff(1:num_rois_2, match_2);

else
    ROIsIDs{1} = ROIsIDs{1}(:)'; % making sure we are [1xN]
    ROIsIDs{2} = ROIsIDs{2}(:)'; % making sure we are [1xN]
    HMatch_results.matched_pairs = [ROIsIDs{1}(match_1) ; ROIsIDs{2}(match_2)];
    HMatch_results.nonmatched_ROIs{1,1} = ROIsIDs{1}(setdiff(1:num_rois_1, match_1));
    HMatch_results.nonmatched_ROIs{1,2} = ROIsIDs{2}(setdiff(1:num_rois_2, match_2));

end



end


