%% Load the greens function produced from Pylith
%
%   [G,patchcoords] = load_pylith_greensfn(folder,faultfilename,Gfilename)
%   - G = [G_horizontal; G_vertical]
%
%   Last modified on 24-Sep-2025
%   by Jeng Hann, Chong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G,lockedpatchcoords,dist,c_disp_ft_sort,c_slip_ft_sort,c_trac_ft_sort] = load_pylith_greensfn(folder,faultfilename,Gfilename)
 

    %%% Load input %%%
    hdf5ft = [folder,faultfilename]; % read the geometry, displacement and greens function for the locked fault 
    
    %%% Read the backslip and sig0 displacement (assuming 1 meter slip) %%%
    hdf5greens = [folder,Gfilename]; % the backslip's greens function and geometry of the rest of the slab top
    

    % For locked fault
    cellsft = h5read(hdf5ft, '/viz/topology/cells');
    verticesft = h5read(hdf5ft, '/geometry/vertices');
    c_slip_ft = h5read(hdf5ft,'/vertex_fields/slip');
    c_trac_ft = h5read(hdf5ft,'/vertex_fields/traction_change');
        
    % For green's function
    cellsft_grns = h5read(hdf5greens, '/viz/topology/cells');
    verticesft_grns = h5read(hdf5greens, '/geometry/vertices');
    c_disp_ft = h5read(hdf5greens,'/vertex_fields/displacement'); % greens function output (assumed we used 1 meter slip to produce this displacement)
    
    %%% Starting the variables organization %%%
    % Make sure locked depth is only from LOCKED FAULT PATCH 
    % if you read the vertices directly from the locked patch, you don't need to change anything
    lockeddepth = min(verticesft(2,:))+0;
    
    % Sorting of stress kernel based on fault depth
    [sorty,I_sortreceiver] = sort(verticesft(2,:),2,'descend'); % need to sort them based on depth
    sortx = verticesft(1,I_sortreceiver); % sort the x-axis
    
    % Portion with greens function (30 kms default, change accordingly)
    locked_I_sortreceiver = I_sortreceiver(sorty>=lockeddepth); % get only the index from the locked portion with greens function (30 kms default, change accordingly)
    
    % sort only the 'locked' vertices (used as impulses) based on locked depth
    % (change locking depth accordingly)
    [sort_cols,I_sortslip] = sort(verticesft(2,verticesft(2,:)>=lockeddepth),2,'descend');
    
    % Sorted position of each patches for locked fault
    sortx_locked = sortx(1:length(locked_I_sortreceiver));
    sorty_locked = sorty(1:length(locked_I_sortreceiver));
    vx = sortx_locked(1:end-1)'; % removing last slip on fault that is not used
    vy = sorty_locked(1:end-1)'; % removing last slip on fault that is not used
    
    lockedpatchcoords = [vx,vy]; % locked patch coordinates

    %%% For backslip (locked patch) %%%
    % Make into a matrix
    [dist,distI] = sort(verticesft_grns(1,:)/1e3); % sorting the distance order [output in km] using GNSS on surface
    [depth,depthI] = sort(verticesft(2,:)/1e3,'descend'); % sorting the depth order [output in km] using locked fault patches
    
    % Final sorting uses both sort indices
    c_disp_ft_sort      = c_disp_ft(:,distI,depthI); % sort the green's function due to slip impose on the locked fault
    c_slip_ft_sort      = c_slip_ft(:,I_sortreceiver,I_sortslip);
    c_trac_ft_sort      = c_trac_ft(:,I_sortreceiver,I_sortslip);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PART II: Running Greens function from Pylith
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We want to get sliprate with 810x1 in the end; but we have both
    % horizontal and vertical displacement, so we vectorize the two displacement
    
    % G matrix for stacked locked horizontal and vertical displacements
    % (removing the last points here)
    G_h = squeeze(c_disp_ft_sort(1,:,1:end-1));  % size: 349 × 811
    G_v = squeeze(c_disp_ft_sort(2,:,1:end-1));  % size: 349 × 811
    
    G = [G_h; G_v];  % size: 698 × 811 % make the horizontal and vertical together [hori up, vert below]
    


end