function [roi_XY_ellipse] = tranform_ROI_2ellipse(roi_XY, size_constrain)

% roi_XY: [N x 2] where N is the number of points of shape


% possible REF: 
% http://www.hep.princeton.edu/mumu/target/Yan/ellipse_fit.pdf
% https://www.cs.cornell.edu/cv/OtherPdf/Ellipse.pdf 


N = size(roi_XY, 1);

% Kind of center data ....
center_init = sum(roi_XY, 1)/size(roi_XY,1);

x = roi_XY(:,1) - center_init(1);
y = roi_XY(:,2) - center_init(2);


% METHOD : Fit with Least Squares
% ---------------------------------
% 
% Solve: A*x^2 + B*x*y + C*y^2 + D*x + E*y + F = 0
% where we divided the ellipse equation by -F

Z = [x.^2 , x.*y , y.^2 , x, y];

param_fit = num2cell(pinv(Z) * ones(N,1));

[A,B,C,D,E] = deal(param_fit{:});

% Get the geometric parameters

% Tilt (counter clockwise, in radians)
if abs(B) < 1e-6
    alpha0 = 0;
else
    alpha0 = 1/2 * acot( (C-A)/B );
end

cos_alpha_2 = cos(alpha0)^2;
sin_alpha_2 = sin(alpha0)^2;

% Recompute the params without the tilt...
[A,C,D,E] = deal(...
    A*cos_alpha_2 + B*cos(alpha0)*sin(alpha0) + C*sin_alpha_2 ,...
    A*sin_alpha_2 - B*cos(alpha0)*sin(alpha0) + C*cos_alpha_2,...
    D*cos(alpha0) + E*sin(alpha0),...
    -D*sin(alpha0) + E*cos(alpha0) );

% Center of ellipse - non-tilted ellipse
center_ellipse = [-D/(2*A) ; -E/(2*C)];

% Get the a and b axis
W = 1 + E^2/(4*C) + D^2/(4*A);
[a,b] = deal( sqrt( abs(W/A) ), sqrt( abs(W/C) ) );


% Get which is the longest axis
if a>b
    thetaE_axis = 0;
else
    thetaE_axis = pi/2;
end

% Check size constraints - cannot vary more than % from original size

% > Get all vectors between points
Vpts = roi_XY - permute(repmat(roi_XY, [1,1,N]),[3,2,1]);

% > Get the length
Vpts_len = sqrt(sum(Vpts.^2, 2)); % [N x 1 x N]

% > Find the longest distance between points, and use it a max
l_max = max(Vpts_len(:)) /2;

a = min(a, l_max*size_constrain);
b = min(b, l_max*size_constrain);


% Set the rotation matrix 

R = [cos(alpha0) , -sin(alpha0) ; sin(alpha0) , cos(alpha0)];


% Finally, get the coorfinates of the shape transformed into ellipse

theta = linspace(0, 2*pi*(2*N-1)/(2*N), 2*N);

roi_XY_ellipse = R * (center_ellipse(:,ones(2*N,1)) + [a*cos(theta) ; b*sin(theta)]);

roi_XY_ellipse = center_init(ones(2*N,1),:) + roi_XY_ellipse';


% Final sanity check, make sure the ellipse is not too shifted 

% Find the unit vector of the axis of the ellipse
vE_axis_norm = R  * [a*cos(thetaE_axis) ; b*sin(thetaE_axis)] / max(a,b); % [2 x 1]

% Project the original points
Proj_roi_vE_axis = roi_XY * vE_axis_norm;

% Get the min/max
min_proj_roi_vE_axis = min(Proj_roi_vE_axis);
max_proj_roi_vE_axis = max(Proj_roi_vE_axis);

% Project the min/max points of the ellipse
xyE_1 = center_init + (R * (center_ellipse + [a*cos(thetaE_axis) ; b*sin(thetaE_axis)]))';
xyE_2 = center_init + (R * (center_ellipse + [a*cos(thetaE_axis+pi) ; b*sin(thetaE_axis+pi)]))';
Proj_E_vE_axis_1 = xyE_1 * vE_axis_norm;
Proj_E_vE_axis_2 = xyE_2 * vE_axis_norm;
min_projE_vE_axis = min([Proj_E_vE_axis_1, Proj_E_vE_axis_2]);
max_projE_vE_axis = max([Proj_E_vE_axis_1, Proj_E_vE_axis_2]);

% Re-align the ellpse to match
if abs(min_proj_roi_vE_axis-min_projE_vE_axis) > 2*abs(max_proj_roi_vE_axis-max_projE_vE_axis)
    mean_xyshift = (max_proj_roi_vE_axis-max_projE_vE_axis) * vE_axis_norm;
elseif 2*abs(min_proj_roi_vE_axis-min_projE_vE_axis) < abs(max_proj_roi_vE_axis-max_projE_vE_axis)
    mean_xyshift = (min_proj_roi_vE_axis-min_projE_vE_axis) * vE_axis_norm;
else
    mean_xyshift = mean([min_proj_roi_vE_axis-min_projE_vE_axis , max_proj_roi_vE_axis-max_projE_vE_axis]) * vE_axis_norm;
end
roi_XY_ellipse = roi_XY_ellipse + repmat(mean_xyshift', [2*N,1]);

end



% % ----------------------------
% 
% % METHOD 1 HOMEMADE....
% % --------------------
% 
% % Compute centroid
% c = sum(roi_XY, 1)/size(roi_XY,1);
% 
% % Get the vectors
% V = [roi_XY(:,1)-c(1), roi_XY(:,2)-c(2)];
% 
% % Get the length
% V_len = sqrt(sum(V.^2, 2)); % [N x 1]
% 
% % Find the longest axis
% [a, idx_a] = max(V_len);
% 
% % Find the corresponding angle (i.e. rotation for the ellipse)
% % theta0 = mod(atan2(V(idx_a,2), V(idx_a,1)), 2*pi); % USING atan2 (use 4 quandrants)
% theta0 = mod(atan(V(idx_a,2)/V(idx_a,1)), 2*pi); % USING atan
% 
% % Estimate b : look at closest points at +/-90 deg - > find the data points
% % with max dot product
% V_norm = V ./ V_len(:,ones(1,2));
% V_Dot_90 = abs( V_norm * [cos(theta0+pi/2); sin(theta0+pi/2)] ); % [N * 1]
% [~, idxMax_Dot_90] = sort(V_Dot_90, 'descend');
% b = mean(V_len(idxMax_Dot_90(1:2)));
% % [~, idxMax_Dot_plus90] = sort(V./V_len(:,ones(1,2)) * [cos(theta0+pi/2); sin(theta0+pi/2)], 'descend');
% % [~, idxMax_Dot_minus90] = sort(V./V_len(:,ones(1,2)) * [cos(theta0-pi/2); sin(theta0-pi/2)], 'descend');
% % b = max( mean(V_len(idxMax_Dot_plus90(1:2))) , mean(V_len(idxMax_Dot_minus90(1:2))));
% % 
% % 
% % % Set the rotation matrix 
% % 
% % R = [cos(theta0) , -sin(theta0) ; sin(theta0) , cos(theta0)];
% % 
% % % Finally, get the coorfinates of the shape transformed into ellipse
% % 
% % theta = linspace(0, 2*pi*(2*N-1)/(2*N), 2*N);
% % 
% % roi_XY_ellipse = R * [a*cos(theta) ; b*sin(theta)];
% % 
% % roi_XY_ellipse = c(ones(2*N,1),:) + roi_XY_ellipse';
% % 


% ----------------------------

% METHOD 2 -- / TO DOUBLE CHECK
% --------------------
% % Get all vectors between points
% Vpts = roi_XY - permute(repmat(roi_XY, [1,1,N]),[3,2,1]);
% Vpts = permute(Vpts, [2,1,3]);
% 
% % Get the length
% Vpts_len = squeeze(sqrt(sum(Vpts.^2, 1))); % [N x N]
% 
% % Find the longest axis
% [a, idx_a] = max(Vpts_len(:));
% a = a/2;
% % [idx_a1,idx_a2] = ind2sub([N,N],idx_a);
% 
% % Find the corresponding angle (i.e. rotation for the ellipse)
% % CHECK USING atan2 (use 4 quandrants)
% theta0 = mod(atan(Vpts(2,idx_a)/Vpts(1,idx_a)), 2*pi); % USING atan
% 
% % Compute centroid cordinates
% % c = [(roi_XY(idx_a1,1)+roi_XY(idx_a2,1))/2, (roi_XY(idx_a1,2)+roi_XY(idx_a2,2))/2];
% c = sum(roi_XY, 1)/size(roi_XY,1);
% 
% % Get all vectors with centroid
% Vc = [roi_XY(:,1)-c(1), roi_XY(:,2)-c(2)];
% 
% % Estimate b : look at closest points at +/-90 deg
% Vc_len = sqrt(sum(Vc.^2, 2));
% Vc_norm = Vc ./ Vc_len(:,ones(1,2));
% Vc_Dot_90 = abs( Vc_norm * [cos(theta0+pi/2); sin(theta0+pi/2)] ); % [N * 1]
% [~, idxMax_Dot_90] = sort(Vc_Dot_90, 'descend');
% b = mean(Vc_len(idxMax_Dot_90(1:2)));
% 
% 
% % Set the rotation matrix 
% 
% R = [cos(theta0) , -sin(theta0) ; sin(theta0) , cos(theta0)];
% 
% % Finally, get the coorfinates of the shape transformed into ellipse
% 
% theta = linspace(0, 2*pi*(2*N-1)/(2*N), 2*N);
% 
% roi_XY_ellipse = R * [a*cos(theta) ; b*sin(theta)];
% 
% roi_XY_ellipse = c(ones(2*N,1),:) + roi_XY_ellipse';



