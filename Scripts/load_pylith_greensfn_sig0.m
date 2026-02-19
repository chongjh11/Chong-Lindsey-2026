%% Load the greens function produced from Pylith
%
%   [G,patchcoords] = load_pylith_greensfn(folder,faultfilename,Gfilename)
%   - G = [G_horizontal; G_vertical]
%
%   Last modified on 24-Sep-2025
%   by Jeng Hann, Chong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G_sig0,lockedpatchcoords,dist,c_disp_ft_sort_sig0,c_slip_ft_sort_sig0,c_trac_ft_sort_sig0] = load_pylith_greensfn_sig0(folder_sig0,faultfilename,Gfilename)
 
    %%% Load input %%%
    % Reading slabtop (all)
    hdf5ftl_sig0 = [folder_sig0,faultfilename];
    
    %%% Read the backslip and sig0 displacement (assuming 1 meter slip) %%%
    hdf5greens_sig0 = [folder_sig0,Gfilename]; % the backslip's greens function for surface displacement
    
    % For slabtop all
    c_slip_ft_sig0 = h5read(hdf5ftl_sig0,'/vertex_fields/slip');
    c_trac_ft_sig0 = h5read(hdf5ftl_sig0,'/vertex_fields/traction_change');
    cellsft_sig0 = h5read(hdf5ftl_sig0, '/viz/topology/cells');
    verticesft_sig0 = h5read(hdf5ftl_sig0, '/geometry/vertices');

    % For green's function (of ESPM, not locked fault)
    cellsft_sig0_disp = h5read(hdf5greens_sig0, '/viz/topology/cells');
    verticesft_sig0_disp = h5read(hdf5greens_sig0, '/geometry/vertices');
    c_disp_ft_sig0 = h5read(hdf5greens_sig0,'/vertex_fields/displacement'); % greens function (assumed we used 1 meter slip to produce this displacement)

    
    %%% Starting the variables organization %%%
    % Make sure locked depth is only from LOCKED FAULT PATCH 
    % if you read the vertices directly from the locked patch, you don't need to change anything
    lockeddepth = min(verticesft_sig0(2,:))+0;
    
    % Sorting of stress kernel based on fault depth
    [sorty,I_sortreceiver_sig0] = sort(verticesft_sig0(2,:),2,'descend'); % need to sort them based on depth
    sortx = verticesft_sig0(1,I_sortreceiver_sig0); % sort the x-axis
    
    % Portion with greens function (30 kms default, change accordingly)
    locked_I_sortreceiver = I_sortreceiver_sig0(sorty>=lockeddepth); % get only the index from the locked portion with greens function (30 kms default, change accordingly)
    
    % sort only the 'locked' vertices (used as impulses) based on locked depth
    % (change locking depth accordingly)
    [sort_cols,I_sortslip] = sort(verticesft_sig0(2,verticesft_sig0(2,:)>=lockeddepth),2,'descend');
    
    % Sorted position of each patches for locked fault
    sortx_locked = sortx(1:length(locked_I_sortreceiver));
    sorty_locked = sorty(1:length(locked_I_sortreceiver));
    vx = sortx_locked(1:end-1)'; % removing last slip on fault that is not used
    vy = sorty_locked(1:end-1)'; % removing last slip on fault that is not used
    
    
    %%% For backslip (locked patch) %%%
    % Make into a matrix
    [dist,distI] = sort(verticesft_sig0_disp(1,:)/1e3); % sorting the distance order [output in km] using GNSS on surface
    % [depth,depthI] = sort(verticesft_sig0(2,:)/1e3,'descend'); % sorting the depth order [output in km] using locked fault patches
    
    % Final sorting uses both sort indices
    c_disp_ft_sort_sig0 = c_disp_ft_sig0(:,distI); % sig0 on locked fault in ESPM
    c_slip_ft_sort_sig0 = c_slip_ft_sig0(:,I_sortreceiver_sig0); % sig0 on locked fault due to ESPM
    c_trac_ft_sort_sig0 = c_trac_ft_sig0(:,I_sortreceiver_sig0); % sig0 on locked fault due to ESPM

    lockedpatchcoords = [vx,vy];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PART II: Running Greens function from Pylith
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We want to get sliprate with 810x1 in the end; but we have both
    % horizontal and vertical displacement, so we vectorize the two displacement
    
    % G matrix for stacked locked horizontal and vertical displacements
    % (removing the last points here)
    G_h_sig0 = c_disp_ft_sort_sig0(1,:)';  % size: 349 × 811
    G_v_sig0 = c_disp_ft_sort_sig0(2,:)';  % size: 349 × 811
    
    G_sig0 = [G_h_sig0; G_v_sig0]; % G for sig0 background relate to surface
    


end