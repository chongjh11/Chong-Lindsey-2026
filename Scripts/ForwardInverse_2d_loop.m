%%  Looping to get the forward and inverse model together from a model from Pylith
%   This script creates the fake data (from forward models) and do the inversion from the already
%   produced green's function (from inverse models) and then compare the
%   slip on the fault.
%
%   Inversion only focuses on backslip (didn't use sig0) other than for
%   plotting.
%   
%   NOTES:
%   - You first need the 2 displacements, slip rates, and stressing rates 
%   from the Pylith forward modeling outputs
%   - Stopped on August 2025 because there was some problems with the
%   results (doesn't fully recover input) but the code works (we believe)
%   - You can change the greens function and deformation/slip rate from
%   another model to test the biasness
%
%   TIPS:
%   - We use this a lot to have a forward modeling from an accretionary
%   prism and use the homogeneous model for inversion
%
%   Last modified on 25-Nov-2025
%   by Jeng Hann, Chong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all

% Put in the path of the forward models (the models you want the displacement from)
%%% CHANGE ACCORDINGLY %%%
allfolds = readtable('./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/10deg40thc/listmatwed10deg.txt', 'ReadVariableNames', false,'Delimiter','.');


for vi = 1:size(allfolds,1) % run each folder for forward modeling; E.g., [1:10]

    for li = 1:10 % how many runs do you want for noises? This is useful for plotting the uncertainty ranges; E.g., [1:10]
                    % this is the number of runs per model
    
    % clear all
    % close all
    clc
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PART I: Forward modeling %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load the paths of Pylith output [CHANGE ACCORDINGLY]
    naming = char(table2cell(allfolds(vi,1))); % for looping the names
    folder_sig0_fwd = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/10deg40thc/elastic_slab_sig0_',naming];
    folder_fwd = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/10deg40thc/elastic_slab_',naming];
    
    % Locking pattern [meters] for input model (forward model) [CHANGE ACCORDINGLY]
    % For 10 deg
    lockstart   = -7.3e3; % below is the fully locked depth (up dip of locked patch) (Dc)
    lockend     = -12e3; % above is the fully locked depth (down dip of locked patch) (Dt)
    driven      = -16e3; % driven creep, slip below this depth (this is downdip of a transition locking depth) (Dl)
    % For 20 deg
    % lockstart   = -15e3; % below is the fully locked depth (up dip of locked patch) (Dc)
    % lockend     = -25e3; % above is the fully locked depth (down dip of locked patch) (Dt)
    % driven      = -34e3; % driven creep, slip below this depth (this is downdip of a transition locking depth) (Dl)
    Vpl         = 0.04; % m/yr (incoming plate rate/far field rate)


    % Don't need to change these below unless needed
    faultfilename_fwd = '/step_slab_greensfns-fault_slabtop_lock.h5'; % load the files
    faultfilename_sig0_fwd = '/step_slab_greensfns_sig0-fault_slabtop.h5'; % load the files
    Gfilename_fwd = '/step_slab_greensfns-boundary.h5'; % load the files
    Gfilename_sig0_fwd = '/step_slab_greensfns_sig0-boundary.h5'; % load the files
    
    % Load G for forward model (don't need to change)
    [G_fwd,faultcoords_fwd,dist_fwd,c_disp_ft_sort_fwd,c_slip_ft_sort_fwd,c_trac_ft_sort_fwd] = load_pylith_greensfn(folder_fwd,faultfilename_fwd,Gfilename_fwd);
    [G_sig0_fwd,faultcoords_sig0_fwd,dist_sig0_fwd,c_disp_ft_sort_sig0_fwd,c_slip_ft_sort_sig0_fwd,c_trac_ft_sort_sig0_fwd] = load_pylith_greensfn_sig0(folder_sig0_fwd,faultfilename_sig0_fwd,Gfilename_sig0_fwd);
    
    % Load the locking pattern parameters (need stress constrain as input)
    [input_sliprate_back,input_sliprate_ESPM,~] = get_sliprate_stressconstrain(lockstart, lockend, driven, Vpl, c_trac_ft_sort_fwd, c_trac_ft_sort_sig0_fwd, faultcoords_fwd);
    
    % Get the displacement if sliprate is from backslip
    % disp_vert = (squeeze(c_disp_ft_sort_fwd(2,:,1:end-1)))*Vpl*(ones(size(input_sliprate_ESPM))); % col is receiver, row is input slip
    % disp_hori = (squeeze(c_disp_ft_sort_fwd(1,:,1:end-1)))*Vpl*(ones(size(input_sliprate_ESPM))); % col is receiver, row is input slip
    
        % d_combined_interseismic_fwd = G_sig0_fwd*Vpl + G_fwd*input_sliprate_ESPM; % Get the observed data (or input)
        % Gvert = vertexftgs(2,:,1:end-1);
        % vertG = reshape(Gvert,size(Gvert,2),size(Gvert,3));
        % inputdispvert_fwd = vertG*input_sliprate_ESPM; % reshape this to 100x103
        d_combined_interseismic_fwd = G_sig0_fwd*Vpl + G_fwd*input_sliprate_ESPM; % Get the observed data (or input)
    
        inputdisphori_fwd = d_combined_interseismic_fwd(1:length(d_combined_interseismic_fwd)/2);
        inputdispvert_fwd = d_combined_interseismic_fwd((length(d_combined_interseismic_fwd)/2)+1:end);
    
    % Get additional information [OPTIONAL]
    KK_fwd = squeeze(c_trac_ft_sort_fwd(2,1:end-1,1:end-1));
    KK_sig0_fwd = squeeze(c_trac_ft_sort_sig0_fwd(2,1:size(KK_fwd,1)));
    
    
    % Plot some surface displacement
    figure(1),clf,hold on
        plot(G_fwd*input_sliprate_ESPM,'kx') % interseismic (includes backslip and plate rate)
        plot(G_fwd*(Vpl - input_sliprate_ESPM)) % corrected to get backslip 
    
        grid on
        legend('Interseismic (has plate rate)','Backslip (due to fault only)')
        xlabel('Patches x2')
        ylabel('Displacement [m]')
        set(gca,'fontsize',15)
    
    
        % % For displacement on surface - used for different geometries
        % hdf5gs = [folder_fwd,'/step_slab_greensfns-gnss_stations.h5'];
        % 
        % cellsgs = h5read(hdf5gs, '/viz/topology/cells');
        % verticesgs = h5read(hdf5gs, '/geometry/vertices');
        % vertexftgs = h5read(hdf5gs, '/vertex_fields/displacement');
        % 
        % 
        % for aa = 1:size(vertexftgs,3)
        %     all_vert(aa,:) = vertexftgs(2,:,aa);
        %     all_hori(aa,:) = vertexftgs(1,:,aa);
        % end
        % sum_vert = sum(all_vert);
        % sum_hori = sum(all_hori);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PART II: Inverse modeling %%%
    % [CHANGE ACCORDINGLY]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load the paths of Pylith output (get greens function from) [CHANGE ACCORDINGLY]
    % A common example we use is the homogeneous ESPM's Greens function
    folder_sig0_inv = './Deformation/Models_final_slab_thickness_2d/elastic_slab_sig0_mesh_tri_slab_10deg_40thc_final_v1'; % CHANGE ACCORDINGLY
    folder_inv = './Deformation/Models_final_slab_thickness_2d/elastic_slab_mesh_tri_slab_10deg_40thc_final_v1'; % CHANGE ACCORDINGLY
    % folder_sig0_inv = folder_sig0_fwd; % if you want the same green's function
    % folder_inv = folder_fwd;

    % Don't need to change these below unless needed
    faultfilename_inv = '/step_slab_greensfns-fault_slabtop_lock.h5';
    faultfilename_sig0_inv = '/step_slab_greensfns_sig0-fault_slabtop.h5';
    Gfilename_inv = '/step_slab_greensfns-boundary.h5';
    Gfilename_sig0_inv = '/step_slab_greensfns_sig0-boundary.h5';
    
    % Load G for inverse model
    [G_inv,faultcoords_inv,dist_inv,c_disp_ft_sort_inv,c_slip_ft_sort_inv,c_trac_ft_sort_inv] = load_pylith_greensfn(folder_inv,faultfilename_inv,Gfilename_inv);
    [G_sig0_inv,faultcoords_sig0_inv,dist_sig0_inv,c_disp_ft_sort_sig0_inv,c_slip_ft_sort_sig0_inv,c_trac_ft_sort_sig0_inv] = load_pylith_greensfn_sig0(folder_sig0_inv,faultfilename_sig0_inv,Gfilename_sig0_inv);
    
    % Load the locking pattern parameters (need stress constrain as input)
    [inv_sliprate_back,inv_sliprate_ESPM,~] = get_sliprate_stressconstrain(lockstart, lockend, driven, Vpl, c_trac_ft_sort_inv, c_trac_ft_sort_sig0_inv, faultcoords_inv);
    
    
      % % For displacement on surface - used for different geometries
        % hdf5gs = [folder_inv,'/step_slab_greensfns-gnss_stations.h5'];
        % 
        % cellsgs = h5read(hdf5gs, '/viz/topology/cells');
        % verticesgs = h5read(hdf5gs, '/geometry/vertices');
        % vertexftgs = h5read(hdf5gs, '/vertex_fields/displacement');
        % 
        % 
        % for aa = 1:size(vertexftgs,3)
        %     all_vert(aa,:) = vertexftgs(2,:,aa);
        %     all_hori(aa,:) = vertexftgs(1,:,aa);
        % end
        % sum_vert_inv = sum(all_vert);
        % sum_hori_inv = sum(all_hori);
    
    
    % get unique input locations
    [dist_fwd_unique,I_uniq_in,Iuniq_out] = unique(dist_fwd);
    
    % interpolate forward model output to inverse modeling locations
    inputdisphori_interp = interp1(dist_fwd_unique',inputdisphori_fwd(I_uniq_in),dist_inv');
    inputdispvert_interp = interp1(dist_fwd_unique',inputdispvert_fwd(I_uniq_in),dist_inv');
    
    % combine horiz and vert
    d_combined_interseismic_interp = [inputdisphori_interp; inputdispvert_interp];
    
    % Add gaussian noise to synthetic input
    noise_amp = 0.001; % change accordingly
    d_combined_interseismic_inv = d_combined_interseismic_interp + noise_amp*randn(size(d_combined_interseismic_interp));
    
    % Remove the deep slip component from the displacement so that we just have shallow slip left
    % Use the displacement from the inverse model, Vpl because of the long term plate rate
        Vpl_vec = Vpl * ones(size(inv_sliprate_ESPM)); % set up the plate velocity for locked patches
        % add the backslip correction using G_inv and G_sig0_inv (because we don't know G_fwd)
        d_no_locking = G_sig0_inv*Vpl + G_inv*Vpl_vec; % displacement when there is full creep to surface
        d_combined_interseismic_inv =  d_no_locking - d_combined_interseismic_inv; % get backslip here by removing the interseismic of homogeneous (inverse) by the forward modeling
        
        inputdisphori_inv = d_combined_interseismic_inv(1:length(d_combined_interseismic_inv)/2);
        inputdispvert_inv = d_combined_interseismic_inv((length(d_combined_interseismic_inv)/2)+1:end);
    
        % figure;plot(G_sig0_inv*Vpl); hold on; plot(d_combined,'k');
    
    %%% Do inversion to get the slip rate on the fault - backslash or lsqlin
    KK_inv = squeeze(c_trac_ft_sort_inv(2,1:end-1,1:end-1));
    KK_sig0_inv = squeeze(c_trac_ft_sort_sig0_inv(2,1:size(KK_inv,1)));
    
    
    %%% Regularization with backslash to get slip (L2 norm) %%%
    lambda = 5e-8; % (CHANGE ACCORDINGLY) increase regularization strength for smoothing 
    
    % Sliprate calculation from our given displacement input and G matrix,
    % smoothing included here
    KK_final = lambda*KK_inv;
    pred_ESPM_slip_1 = [G_inv; KK_final]\[d_combined_interseismic_inv; zeros(size(G_inv,2),1)]; % using stress kernel as smoothing
    
    
    %%% Constrained regularization using lsqlin to get slip 
    % [OPTIONAL] Add in parameter for locking depth for KK; deeper than the given locked depth will be "0"
    % Comment if not using this section
    KK_crop = KK_inv;
    stressmaxdepth = -38e3; % [negative for depth in meters] CHANGE ACCORDINGLY
    Iunlocked = (stressmaxdepth>=faultcoords_fwd(:,2)); % select anything below this depth
    KK_crop(Iunlocked,:) = []; % make all those below the locked depth to zero
    
    % Set the upper and lower bounds here, if not just do = "[]"
    lb = zeros(length(KK_sig0_inv),1); % zero slip [meters]
    ub = Vpl + zeros(length(KK_sig0_inv),1); % plate rate [meters]
    
    % Get the stresses only above the locking depth
    KK_sig0_crop = KK_sig0_inv;
    KK_sig0_crop(Iunlocked)=[];
    
    %%% Change your data input restrictions %%%
    % E.g., If you want to invert with horizontal data onshore only, like in
    % real-life or full data onshore and offshore
    %%% (1) fit to all vertical and horizontal onshore and offshore, no stress constraints  %%%
    % fulloffonshore
    % Iland = true(1,length(dist_fwd));
    % pred_ESPM_slip_2 = lsqlin([G_inv; 1*lambda*KK_inv], [d_combined_interseismic_inv; zeros(size(G_inv,2),1)],[],[],[],[],lb,ub); % getting the ESPM slip
    
    %%% (2) fit only to horizontals onshore and offshore, no stress constraints %%%
    % horioffonshore
    % Iland = true(1,length(dist_fwd));
    % pred_ESPM_slip_2 = lsqlin([G_inv(1:end/2,:); 1*lambda*KK_inv], [d_combined_interseismic_inv(1:end/2); zeros(size(G_inv,2),1)],[],[],[],[],lb,ub); % getting the ESPM slip
    % 
    %%% (3) fit to data on full onshore (both vert and hori), no stress constraints %%%
    % fullonshore
    Iland = (dist_fwd>80 & dist_fwd<400); % get x-axis onshore only
    pred_ESPM_slip_2 = lsqlin([G_inv(Iland,:); 1*lambda*KK_inv], [d_combined_interseismic_inv(Iland,:); zeros(size(G_inv,2),1)],[],[],[],[],lb,ub); % getting the ESPM slip

    %%% (4) fit to data on onshore horizontal data only, no stress constraints %%%
    % horionshore
    % Ghori = G_inv(1:end/2,:); % green's function of horizontal only
    % Iland=(dist_fwd>80 & dist_fwd<400);
    % pred_ESPM_slip_2 = lsqlin([Ghori(Iland,:); 1*lambda*KK_inv], [inputdisphori_inv(Iland,:); zeros(size(G_inv,2),1)],[],[],[],[],lb,ub); % getting the ESPM slip
    
    %%% (5) lsqlin constrained slip rate, stress constrained %%%
    % pred_ESPM_slip_2 = lsqlin([G_inv; 1*lambda*KK_inv], [d_combined_interseismic_inv; zeros(size(G_inv,2),1)],KK_crop,KK_sig0_crop,[],[],lb,ub); % getting the ESPM slip
    % pred_ESPM_slip_2 = lsqlin([G_inv; 1*lambda*KK_inv], [d_combined_interseismic_inv; zeros(size(G_inv,2),1)],KK_crop,-KK_sig0_crop*Vpl,[],[],lb,ub); % getting the ESPM slip
    
    %%% (6) fit to data onland only with stress constraints %%%
    % Iland=(dist_fwd>80 & dist_fwd<400);
    % pred_ESPM_slip_2 = lsqlin([G_inv(Iland,:); 1*lambda*KK_inv], [d_combined_interseismic_inv(Iland,:); zeros(size(G_inv,2),1)],KK_crop,-KK_sig0_crop*Vpl,[],[],lb,ub); % getting the ESPM slip
    
    %%% (7) fit to data on land only, with a penalty parameter %%%
    % beta=1e-3;
    % betafit=1; % 1 or 0 - 1 means "no creep", 0 means "all creep
    % Iland=(dist_fwd>50 & dist_fwd<400);
    % pred_ESPM_slip_2 = lsqlin([G_inv(Iland,:); 1*lambda*KK_inv; 1*beta*ones(size(KK_sig0_inv))], [d_combined_interseismic_inv(Iland,:); zeros(size(G_inv,2),1); betafit],[],[],[],[],lb,ub); % getting the ESPM slip
    
    
    % Now using the predicted slip and greens function to get the displacement
    final_pred_disp_ESPM = d_no_locking - G_inv*pred_ESPM_slip_2;
    final_pred_disp_ESPM_hori = final_pred_disp_ESPM(1:length(d_combined_interseismic_inv)/2); % get the horizontal
    final_pred_disp_ESPM_vert = final_pred_disp_ESPM(1+(length(d_combined_interseismic_inv)/2):end); % get the vertical
    
    % For plotting accorindgly
    vx_fwd = faultcoords_fwd(:,1); % number of patches with depth 
    vy_fwd = faultcoords_fwd(:,2); % number of patches with depth 
    
    vx_inv = faultcoords_inv(:,1); % number of patches with depth 
    vy_inv = faultcoords_inv(:,2); % number of patches with depth 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART III: Plot the result
    % In the slip, we subtract Vpl to flip it
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(2); clf, hold on
        plot(-vy_fwd/1e3,Vpl - (Vpl - input_sliprate_ESPM) )
        %plot(-vy,Vpl - pred_ESPM_slip_1)
        plot(-vy_inv/1e3,Vpl - pred_ESPM_slip_2)
        
        grid on
    
        title('Slip rate comparison')
        legend('Input','Lsqlin')
        xlabel('Depth [km]')
        ylabel('Slip rate [m/yr]')
        set(gca,'fontsize',16)
    
    
    %%% Plot deformation and slip rate from ESPM ONLY using stressing kernel here %%%
    f3 = figure(3); %clf; 
    f3.Position = [1321         377         600         600];
    [~, lastPart] = fileparts(folder_fwd);
    escapedTitle = strrep(lastPart, '_', '__');  
    title(escapedTitle)
    
    subplot(3,1,3); hold on % plot the slip rate on the fault
        plot(vx_fwd/1e3,input_sliprate_ESPM,'-k','linewidth',3);
        plot(vx_inv/1e3,Vpl - pred_ESPM_slip_2,'r','linewidth',2)
    
        sgtitle(escapedTitle)
        % title('Slip rate on fault')
        ylabel({'Slip rate on fault','[m/yr]'}); xlabel('Depth from trench [km]')
        legend('Input','Predicted','Location','northwest')
        grid on
        grid minor
        set(gca,'fontsize',16)
    
    subplot(3,1,2); hold on % plot the horizontal due to the slip
        plot(dist_fwd,inputdisphori_fwd,'color',[0.6,0.6,0.6],'linewidth',3) % observed horizontal
        plot(dist_fwd(Iland),inputdisphori_fwd(Iland),'k','linewidth',3) % observed horizontal
        plot(dist_inv,final_pred_disp_ESPM_hori,'r','linewidth',2) % predicted horizontal
    
        % title('Horizontal Deformation')
        ylabel({'Horizontal Deformation','[m/yr]'}); xlabel('Distance from trench [km]')
        legend('Observed data','Land data','Predicted')
        grid on
        ax.XMinorGrid = 'on';
        xlim([-10 220])
        % xlim([min(vx/1e3) max(vx/1e3)])
        % ylim([-0.0 0.05])
        set(gca,'fontsize',16)
    
    s1 = subplot(3,1,1); hold on % plot the vertical due to the slip
        plot(dist_fwd,inputdispvert_fwd,'color',[0.6,0.6,0.6],'linewidth',3) % observed vertical
        plot(dist_fwd(Iland),inputdispvert_fwd(Iland),'k','linewidth',3) % observed vertical
        plot(dist_inv,final_pred_disp_ESPM_vert,'r','linewidth',2) % predicted vertical
    
        % title('Vertical Deformation')
        ylabel({'Vertical Deformation','[m/yr]'}); xlabel('Distance from trench [km]')
        legend('Observed data','Land data','Predicted')
        grid on
        % grid minor
        ax.XMinorGrid = 'on';
        xlim([-10 220])
        % xlim([min(vx/1e3) max(vx/1e3)])
        ylim([-0.015 0.02])
        set(gca,'fontsize',16)
    
    
    %%% Plot the Stressing rates %%%
    f4 = figure(4); clf, hold on
    
        % Observed stresses (from Pylith)
        plot(vy_inv/1e3,KK_inv*pred_ESPM_slip_2,'k-','LineWidth',4)
        
        % Predicted stresses (stress = stress kernel * slip)
        plot(vy_inv/1e3,KK_inv*pred_ESPM_slip_2,'b-','LineWidth',1)
        plot(vy_inv/1e3,KK_sig0_inv*Vpl,'r-','LineWidth',1)
        
        set(gca,'xdir','rev')
        grid on
        title('Stressing rates',escapedTitle)
        xlabel('Depth [km]'); ylabel('Stressing rate')
        legend('Forward model (Pylith) - ESPM','Model fit - ESPM','location','southwest')
        set(gca,'fontsize',15)
    
    
        % Saving the min and max range
        compiled_final_pred_vert(:,li) = final_pred_disp_ESPM_vert;
        compiled_final_pred_hori(:,li) = final_pred_disp_ESPM_hori;
        compiled_final_slip(:,li)      = Vpl - pred_ESPM_slip_2;
        compiled_input_vert(:,li)      = inputdispvert_fwd;
        compiled_input_hori(:,li)      = inputdisphori_fwd;
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% [OPTIONAL] Extract and plot the min and max within each x-axis,
    %%% useful for plotting ranges of min max across transect
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear range_final_slip range_final_pred_vert range_final_pred_hori range_input_pred_vert range_input_pred_hori
    
    % Get the range for horizontal and vertical deformation
    for ci = 1:size(compiled_input_vert,1)
        range_input_pred_vert(ci,:) = [min(compiled_input_vert(ci,:)),max(compiled_input_vert(ci,:))];
        range_input_pred_hori(ci,:) = [min(compiled_input_hori(ci,:)),max(compiled_input_hori(ci,:))];
    
        if ci == size(compiled_input_vert,1)
            range_input_pred_vert(:,size(range_input_pred_vert,2)+1) = dist_fwd';
            range_input_pred_hori(:,size(range_input_pred_vert,2)+1) = dist_fwd';
        end
    end
    
    for ci = 1:size(compiled_final_pred_vert,1)
        range_final_pred_vert(ci,:) = [min(compiled_final_pred_vert(ci,:)),max(compiled_final_pred_vert(ci,:))];
        range_final_pred_hori(ci,:) = [min(compiled_final_pred_hori(ci,:)),max(compiled_final_pred_hori(ci,:))];
        
        if ci == size(compiled_final_pred_vert,1)
            range_final_pred_vert(:,size(range_final_pred_vert,2)+1) = dist_inv';
            range_final_pred_hori(:,size(range_final_pred_vert,2)+1) = dist_inv';
        end
    end
    
    % Get the range for slip on fault
    for cii = 1:size(compiled_final_slip,1)
        range_final_slip(cii,:) = [min(compiled_final_slip(cii,:)),max(compiled_final_slip(cii,:))];
        mean_final_slip(cii,:) = mean(compiled_final_slip(cii,:)); % get the mean for plotting later on
    
        if cii == size(compiled_final_slip,1)
            range_final_slip(:,size(range_final_slip,2)+1) = -vy_inv; % change to "-vy_inv" if you want to plot with depth or "vx_inv" for distance
            range_final_slip(:,size(range_final_slip,2)+1) = vx_inv; % change to "-vy_inv" if you want to plot with depth or "vx_inv" for distance
        end
    end
    
    % Redo the variable for plotting the "fill" for uncertainty range    
    XXvert_inp = [range_input_pred_vert(:,end);flipud(range_input_pred_vert(:,end))];
    YYvert_inp = [range_input_pred_vert(:,1);flipud(range_input_pred_vert(:,2))];
    XXhori_inp = [range_input_pred_hori(:,end);flipud(range_input_pred_hori(:,end))];
    YYhori_inp = [range_input_pred_hori(:,1);flipud(range_input_pred_hori(:,2))];
    XXvert = [range_final_pred_vert(:,end);flipud(range_final_pred_vert(:,end))];
    YYvert = [range_final_pred_vert(:,1);flipud(range_final_pred_vert(:,2))];
    XXhori = [range_final_pred_hori(:,end);flipud(range_final_pred_hori(:,end))];
    YYhori = [range_final_pred_hori(:,1);flipud(range_final_pred_hori(:,2))];
    XXslip = [range_final_slip(:,end);flipud(range_final_slip(:,end))];
    YYslip = [range_final_slip(:,1);flipud(range_final_slip(:,2))];

    % Plot the deformation of the uncertainty range (if only run, it won't plot)
    f5 = figure(5);clf
    f5.Position = [1     377   600   600];

        subplot(3,1,1); hold on % plot vertical 
        % h11 = fill(XXvert_inp,YYvert_inp,[0.1 0.1 0.1], 'EdgeColor','none'); % input w/range
        plot(dist_fwd,inputdispvert_fwd,'color',[0.6,0.6,0.6],'linewidth',3) % observed vertical
        h11 = plot(XXvert_inp(Iland),YYvert_inp(Iland),'k','linewidth',3); % input vertical
        h1 = fill(XXvert,YYvert,[0.8 0.1 0.1], 'EdgeColor','none');
    
        xlim([-10 220]);ylim([-0.015 0.015])
        set(gca,'fontsize',15)
        % set(h11,'FaceAlpha',0.4)                        
        set(h1,'FaceAlpha',0.3)                        
        grid on
        grid minor
        legend('Input','Predicted')
        ylabel({'Vertical Deformation','[m/yr]'}); xlabel('Distance from trench [km]')
    
    
        subplot(3,1,2); hold on % plot horizontal
        % h21 = fill(XXhori_inp,YYhori_inp,[0.1 0.1 0.1], 'EdgeColor','none'); % input w/range
        plot(dist_fwd,inputdisphori_fwd,'color',[0.6,0.6,0.6],'linewidth',3) % observed horizontal
        h21 = plot(XXhori_inp(Iland),YYhori_inp(Iland),'k','linewidth',3); % input horizontal
        h2 = fill(XXhori,YYhori,[0.8 0.1 0.1], 'EdgeColor','none');
    
        xlim([-10 220])
        set(gca,'fontsize',15)
        % set(h21,'FaceAlpha',0.4)                        
        set(h2,'FaceAlpha',0.3)                        
        grid on
        grid minor
        legend('Input','Predicted')
        
        ylabel({'Horizontal Deformation','[m/yr]'}); xlabel('Distance from trench [km]')
        ylim([-0.0 0.045])
        % ylim([-0.0 0.05])
    
    
        subplot(3,1,3); hold on % plot slip rate
        h33 = plot(vx_fwd/1e3,input_sliprate_ESPM,'-k','linewidth',3); % input
        % h33 = plot(-vy_fwd/1e3,input_sliprate_ESPM,'-k','linewidth',3); % input
        h3 = fill(XXslip/1e3,YYslip,[0.8 0.1 0.1], 'EdgeColor','none'); % predicted (change the x-axis to depth or distance above, not here)
        ylabel({'Slip rate on fault','[m/yr]'}); xlabel('Distance from trench [km]')
    
        set(h3,'FaceAlpha',0.2)                        
        set(gca,'fontsize',15)
        grid on
        grid minor
        legend('Input','Predicted','location','northwest')
        sgtitle('Onshore vertical and horizontal data for inversion only')
        xlim([0 110])

%%% Getting the area above the curve to give a more consistent way of looking at the different model fitting
% Get the area of each curves first

    % Get the difference area of the two curves, convert everything to meters
    f6 = figure(6); clf
    subplot(3,1,1); hold on
        % pltfwd_area = area(-vy_fwd/1e3,(input_sliprate_ESPM-Vpl)); % input
        pltfwd_area = area(vx_fwd/1e3,(input_sliprate_ESPM-Vpl)); % input
        pltinv_area = area(XXslip(1:length(XXslip)/2)/1e3,(mean_final_slip-Vpl)); % using the mean of prediction rather than the range
        
        % Slip calculation    
        % fwd_area = trapz(-vy_fwd,Vpl-input_sliprate_ESPM); % input
        fwd_area = trapz(vx_fwd,Vpl-input_sliprate_ESPM); % input
        inv_area = trapz(XXslip(1:length(XXslip)/2),Vpl-mean_final_slip); % predicted
        inv_area_min = trapz(XXslip(1:length(XXslip)/2),Vpl-range_final_slip(:,1)); % min slip
        inv_area_max = trapz(XXslip(1:length(XXslip)/2),Vpl-range_final_slip(:,2)); % max slip
        
        % Get the difference between area (this is m^2/yr)
        % if negative area it means predicted area larger than input
        diff_area(vi,1) = (inv_area-fwd_area); % meters/year X meters = m^2/yr (Mean)
        diff_area(vi,2) = (inv_area_min-fwd_area); % meters/year X meters = m^2/yr (Min)
        diff_area(vi,3) = (inv_area_max-fwd_area); % meters/year X meters = m^2/yr (Max)
        pct_area(vi,1) = diff_area(vi,1)/fwd_area; % in % mean
        pct_area(vi,2) = diff_area(vi,2)/fwd_area; % in % min 
        pct_area(vi,3) = diff_area(vi,3)/fwd_area; % in % max
    
        ylabel({'Slip rate on fault','[m/yr]'}); xlabel('Distance from trench [km]')
        sgtitle(['InputArea: ',num2str(fwd_area),' PredArea: ',num2str(inv_area),' DiffArea =  ',num2str(diff_area(vi,1))])
        grid on
        set(gca,'fontsize',16)
        legend('Input','Predicted','Location','northwest')
        
        pltinv_area.FaceAlpha = 0.5;
        pltfwd_area.FaceAlpha = 1;

    
    %%% Save figures? %%%
    % For each folder
    savi = 1;
    if savi ~= 0
        newpath = [pwd,'/',strrep(folder_fwd, './', '')];
        fprintf('Saving figure...')
    
        mkdir([newpath, filesep, 'Figures'])
        saveas(f5,fullfile(newpath,'/Figures/','Inversion_Homogeneous10deg40thc_fulllocked_fullonshore.png'),'png');  
        saveas(f5,fullfile(newpath,'/Figures/','Inversion_Homogeneous10deg40thc_fulllocked_fullonshore.svg'),'svg');  
 
        saveas(f6,fullfile(newpath,'/Figures/','AreaDiffInversion_Homogeneous10deg40thc_fulllocked_fullonshore.png'),'png');  
        saveas(f6,fullfile(newpath,'/Figures/','AreaDiffInversion_Homogeneous10deg40thc_fulllocked_fullonshore.svg'),'svg');  

        T_rng = table(XXvert,YYvert,XXhori,YYhori, 'VariableNames', {'XXvert[m]','YYvert[m/yr]','XXhori[m]','YYhori[m/yr]'});
        T_slp = table(XXslip,YYslip, 'VariableNames', {'XXslip[m]','YYslip[m/yr]'});
        writetable(T_rng, fullfile((newpath), 'Range_minmax_disp.txt'));
        writetable(T_slp, fullfile((newpath), 'Range_minmax_slip.txt'));


        fprintf('done\n')
        close all

    end

end

fprintf('done looping\n')
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Save the percent area difference as textfiles 
% Default uses all the names from "allfolds, you might get error otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = table(diff_area(:,1),diff_area(:,2),diff_area(:,3),pct_area(:,1),pct_area(:,2),pct_area(:,3),char(table2cell(allfolds(:,1))), 'VariableNames', {'mean_diff_area[m^2/yr]','min_diff_area[m^2/yr]','max_diff_area[m^2/yr]','mean_%_area[m^2/yr]','min_%_area[m^2/yr]','max_%_area[m^2/yr]','filename'});
% writetable(T, fullfile(fileparts(newpath), 'area_difference_inv_fwd_fulloffonshore.txt'));
writetable(T, fullfile(fileparts(newpath), 'area_difference_inv_fwd_fullonshore.txt'));
% writetable(T, fullfile(fileparts(newpath), 'area_difference_inv_fwd_horioffonshore.txt'));
% writetable(T, fullfile(fileparts(newpath), 'area_difference_inv_fwd_horionshore.txt'));

disp('done saving txt files')