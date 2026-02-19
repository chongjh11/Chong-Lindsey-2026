function [interp_greens] = greensinterpolate(ori_dist,ori_greens,mod_dist)

%%% Greens function interpolation %%%
%
%   NOTES:
%   1. "mod_dist" is the gnss locations
%
%   Last modified on 5-Dec-2025
%   by Jeng Hann Chong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the arrays
nobs_orig = length(ori_dist);
nobs_gnss = length(mod_dist);
nmodel = size(ori_greens,2);

% Fix duplicate zeros on dist_inv (from Pylith, pylith uses two "0" to get the displacement at the fault
dist_inv_fix = ori_dist;
Izero = find(ori_dist == 0);
dist_inv_fix(Izero(1))=-1e-9; % set one of the zeros to negative
dist_inv_fix(Izero(2))=1e-9; % set one of the zeros to positive

interp_greens = zeros(2*nobs_gnss, nmodel); % set up the array for new greens function

for icell = 1:nmodel
    % interpolate this column of ori_greens onto gnss_greens (interp_greens)
    ghoriz = ori_greens(1:end/2, icell);
    gvert  = ori_greens(end/2+1:end, icell);

    ghoriz_gnss = interp1(dist_inv_fix,ghoriz,mod_dist,'spline');
    gvert_gnss  = interp1(dist_inv_fix,gvert, mod_dist,'spline');

    interp_greens(1:end/2, icell) = ghoriz_gnss;
    interp_greens(end/2+1:end, icell) = gvert_gnss;
end

%%% Plot to visualize if needed %%%
% figure
% plot(ori_dist,ori_greens(1:end/2 ,100),'-b.') % original hori
% hold on
% plot(ori_dist,ori_greens(end/2+1:end ,100),'-r.') % original vert
% plot(mod_dist,0*mod_dist,'x') % location of GNSS
% plot(mod_dist,interp_greens(1:end/2 ,100),'--kx') % interp hori
% hold on
% plot(mod_dist,interp_greens(end/2+1:end ,100),'--mx') % interp vert

end