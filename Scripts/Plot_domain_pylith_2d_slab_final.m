%% Plot the domain and GNSS displacement from Green's functions
%
%   NOTES:
%   - This can show us where the each green's function is located and the rates
%   - You can use this to show the order of the green's function too
%   - The traction organization (stress kernel/green's) is [component,all receivers,slip imposed]
%
%   - I_sortreceiver: number of observer patches (includes locked and
%   unlocked)
%   - I_sortslip: number of slip patches (patches that are locked)
%   - This includes several OPTIONAL sections you can plot to show more things
%
%
%   INSTRUCTIONS:
%   1. The directory to run the script should be before the folder you save
%   your different versions of case studies (within MATLAB > Pylith)
%   2. You can run the output from Pylith 
%   3. You can save the textfile version of the displacements (after PART
%   VI) & the figures of displacement and stressses in PART IV
%   4. It produces textfiles that can be used by
%   "Plot_compilation_modelresults.m".
%
%   
%   Last modified on 29-Aug-2025
%   by Jeng Hann, Chong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART 0: Renaming folders
% 
%   Rename the folders if you haven't. You need to rename it based on the
%   list you give, each version has a different property, check before
%   renaming. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs here
% ni = 369:384; % change accordingly
% prifold = 'Models_final_slab_wedge_elastic_2d'; % change accordingly
% secfold = '25GPaWedge_60GPaSlab_30GPaCrust_70GPaMantle'; % change accordingly
% 
% meshnames = readtable('./Deformation/Models_final_slab_wedge_elastic_2d/listmeshwedge_rename.txt','Delimiter','');
% 
% % Renaming here
% for vi = 1:length(ni)
% 
%     oldFolder = fullfile(pwd, './Deformation/', prifold, '/', secfold, strcat('elastic_slab_',string(table2array(meshnames(vi,:)))));
%     newFolder = fullfile(pwd, './Deformation/', prifold, '/', secfold, ['/elastic_slab_v', num2str(ni(vi))]);
% 
%     fprintf('Renaming: "%s" to "%s"',strcat('elastic_slab_',string(table2array(meshnames(vi,:)))),['/elastic_slab_v', num2str(ni(vi))])
%     f=movefile(oldFolder, newFolder); % for elastic_slab
% 
%     oldFolder2 = fullfile(pwd, './Deformation/', prifold, '/', secfold, strcat('elastic_slab_sig0_',string(table2array(meshnames(vi,:)))));
%     newFolder2 = fullfile(pwd,'./Deformation/', prifold, '/', secfold, ['/elastic_slab_sig0_v',num2str(ni(vi))]);
%     fprintf('  Renaming: "%s" to "%s"',strcat('elastic_slab_sig0_',string(table2array(meshnames(vi,:)))),['/elastic_slab_sig0_v',num2str(ni(vi))])
%     f2=movefile(oldFolder2, newFolder2); % for elastic_slab_sig0
% 
%     fprintf('done\n')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART I: Input data
% 
% Normal folder "sl_vX" = backslip (green's function on locked patch)
% Sig0 folder "sig0_sl_vX" = deep creep/background (slip on outside locked)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all

counter = 0; % IGNORE THIS: for saving the textfiles later in the section

% Alternative way of reading many folders
allfolds = readtable('./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/30deg60thc/listmatwed30deg.txt', 'ReadVariableNames', false,'Delimiter','.');
% allfolds = readtable('./Deformation/Models_final_slab_thickness_2d/listmesh.txt', 'ReadVariableNames', false,'Delimiter','.');

for vi = 1:size(allfolds,1)%[1:10] % put in the number versions (e.g., elastic_ longoweslab_vX). This is for the "X" numbers
close all

%%% Setting up %%%

counter = counter + 1; % for saving the textfiles later in the section
 
%%% CHANGE ACCORDINGLY to the input file here after changing the textfile above %%%
% folder_sig0 = ['./Deformation/Models_final_slb_wedge_elastic_2d/1GPaWedge_30GPaSlab_30GPaCrust_70GPaMantle/elastic_slab_sig0_v',num2str(vi)]; 
% folder = ['./Deformation/Models_final_slab_wedge_elastic_2d/1GPaWedge_30GPaSlab_30GPaCrust_70GPaMantle/elastic_slab_v',num2str(vi)];
% folder_sig0 = ['./Deformation/Models_final_slab_elastic_2d/30GPaSlab_30GPaCrust_70GPaMantle/elastic_slab_sig0_v',num2str(vi)]; 
% folder = ['./Deformation/Models_final_slab_elastic_2d/30GPaSlab_30GPaCrust_70GPaMantle/elastic_slab_v',num2str(vi)];
% folder_sig0 = ['./Deformation/Models_final_slab_wedge_elastic_2d/test/elastic_slab_s spuig0_v',num2str(vi)]; 
% folder = ['./Deformation/Models_final_slab_wedge_elastic_2d/test/elastic_slab_v',num2str(vi)];
% folder_sig0 = ['./Deformation/Models_final_slab_elastic_2d/test/elastic_slab_sig0_v',num2str(vi)]; 
% folder = ['./Deformation/Models_final_slab_elastic_2d/test/elastic_slab_v',num2str(vi)];
% folder_sig0 = ['./Deformation/Models_final_slab_matwedge_elastic_2d/10deg40thc/elastic_slab_sig0_v',num2str(vi)]; 
% folder = ['./Deformation/Models_final_slab_matwedge_elastic_2d/10deg40thc/elastic_slab_v',num2str(vi)];
% folder_sig0 = ['./Deformation/Models_final_slab_thickness_2d/elastic_slab_sig0_v',num2str(vi)];
% folder = ['./Deformation/Models_final_slab_thickness_2d/elastic_slab_v',num2str(vi)];
% folder_sig0 = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/20deg40thc/elastic_slab_sig0_v',num2str(vi)];
% folder = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/20deg40thc/elastic_slab_v',num2str(vi)];

naming = char(table2cell(allfolds(vi,1)));
folder_sig0 = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/30deg60thc/elastic_slab_sig0_',naming];
folder = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/30deg60thc/elastic_slab_',naming];
% folder_sig0 = ['./Deformation/Models_final_slab_thickness_2d/elastic_slab_sig0_',naming];
% folder = ['./Deformation/Models_final_slab_thickness_2d/elastic_slab_',naming];

%%% Load input %%%
% hdf5 = [folder,'/step_slab_greensfns-domain.h5'];
hdf5ft = [folder,'/step_slab_greensfns-fault_slabtop_lock.h5'];
% hdf5ft_sig0 = [folder_sig0,'/step_slab_greensfns_sig0-fault_slabtop_lock.h5'];

% Reading slabtop (all)
% hdf5ftl = [folder_sig0,'/step_slab_greensfns_sig0-fault_slabtop_lower.h5'];
hdf5ftl_sig0 = [folder_sig0,'/step_slab_greensfns_sig0-fault_slabtop.h5'];

% Reading slabbot
hdf5sb = [folder_sig0,'/step_slab_greensfns_sig0-fault_slabbot.h5'];
    
    % % For displacement on surface
    % cellsgs = h5read(hdf5gs, '/viz/topology/cells');
    % verticesgs = h5read(hdf5gs, '/geometry/vertices');
    % vertexftgs = h5read(hdf5gs, '/vertex_fields/displacement');

% For locked fault
cellsft = h5read(hdf5ft, '/viz/topology/cells');
verticesft = h5read(hdf5ft, '/geometry/vertices');
c_slip_ft = h5read(hdf5ft,'/vertex_fields/slip');
c_trac_ft = h5read(hdf5ft,'/vertex_fields/traction_change');


% For slabtop all
c_slip_ft_sig0 = h5read(hdf5ftl_sig0,'/vertex_fields/slip');
c_trac_ft_sig0 = h5read(hdf5ftl_sig0,'/vertex_fields/traction_change');
cellsft_sig0 = h5read(hdf5ftl_sig0, '/viz/topology/cells');
verticesft_sig0 = h5read(hdf5ftl_sig0, '/geometry/vertices');

% Make sure locked depth is only from LOCKED FAULT PATCH 
% if you read the vertices directly from the locked patch, you don't need to change anything
lockeddepth = min(verticesft(2,:))+0;


%%% Sorting below here %%%
% Sorting of stress kernel based on fault depth
[sorty,I_sortreceiver] = sort(verticesft(2,:),2,'descend'); % need to sort them based on depth
sortx = verticesft(1,I_sortreceiver); % sort the x-axis

% Portion with greens function (30 kms default, change accordingly)
locked_I_sortreceiver = I_sortreceiver(sorty>=lockeddepth); % get only the index from the locked portion with greens function (30 kms default, change accordingly)

% sort only the 'locked' vertices (used as impulses) based on locked depth
% (change locking depth accordingly)
[sort_cols,I_sortslip] = sort(verticesft(2,verticesft(2,:)>=lockeddepth),2,'descend');

% Sorting based on fault depth
[sorty2,I_sortreceiver_sig0] = sort(verticesft_sig0(2,:),2,'descend'); % need to sort them based on depth
sortx2 = verticesft_sig0(1,I_sortreceiver_sig0); % sort the x-axis

% Portion with greens function (30 kms default, change accordingly)
% locked_I_sortreceiver = I_sortreceiver(sorty>=lockeddepth); % get only the index from the locked portion with greens function (30 kms default, change accordingly)

% sort only the 'locked' vertices (used as impulses) based on locked depth
% (change locking depth accordingly)
% [sort_cols,I_sortslip] = sort(verticesft(2,verticesft(2,:)>=lockeddepth),2,'descend');


% Final sorting uses both sort indices
c_slip_ft_sort = c_slip_ft(:,I_sortreceiver,I_sortslip);
c_trac_ft_sort = c_trac_ft(:,I_sortreceiver,I_sortslip);
c_slip_ft_sort_sig0 = c_slip_ft_sig0(:,I_sortreceiver_sig0); % sig0 on locked fault due to ESPM
c_trac_ft_sort_sig0 = c_trac_ft_sig0(:,I_sortreceiver_sig0); % sig0 on locked fault due to ESPM


% Sorted position of each patches for locked fault
sortx_locked = sortx(1:length(locked_I_sortreceiver));
sorty_locked = sorty(1:length(locked_I_sortreceiver));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PART II: Fake unicycle calculations for the backslip case
% Calculate the stressing rate and slip rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CHANGE HERE %%% - constant backslip: units of meters/yr, scaled by 1/cos(dip)
inputsliprate = 0.04; % [m/yr] change accordingly

% Drop the last vertex because it plots wrongly (Pylith output wrongly maybe due to its calculation at fault edge) 
KK = squeeze(c_trac_ft_sort(2,1:end-1,1:end-1)); % stress kernel on the locked patch
vx = sortx_locked(1:end-1)';
vy = sorty_locked(1:end-1)';

KK_sig0 = squeeze(c_trac_ft_sort_sig0(2,1:size(KK,1)))'; % the sig0 due to ESPM model (influence on locked patch from slabbot etc)

Vpl = inputsliprate*ones(size(vx)); % locked patches slip rates

sig0_backslip = KK*Vpl; % synthetic stress for the backslip case 
sig0_slab = -KK_sig0.*Vpl; % stress calculated from Pylith for ESPM [optional, depends on the input slip in Pylith]


figure(1), hold on
    plot(vx/1e3,sig0_backslip,'k') % just locked patch (one fault model)
    plot(vx/1e3,sig0_slab,'-r') % slab ESPM

    legend('Locked patch only','Includes deep dislocation & slab bottom')
    title('Sig0')
    xlabel('Distance [km]')
    grid on


%%% Define a locking patch and transition locking at depth (CHANGE ACCORDINGLY) %%%
% between "lockend" and "driven" is the locking depth transition

% For 10 deg
% lockstart   = 0e3; % below is the fully locked depth (up dip of locked patch) (Dc)
% lockstart   = -7.3e3; % below is the fully locked depth (up dip of locked patch) (Dc)
% lockend     = -12e3; % above is the fully locked depth (down dip of locked patch) (Dt)
% driven      = -16e3; % driven creep, slip below this depth (this is downdip of a transition locking depth) (Dl)

% For 20 deg
% lockstart   = -0e3; % below is the fully locked depth (up dip of locked patch) (Dc)
lockstart   = -15e3; % below is the fully locked depth (up dip of locked patch) (Dc)
lockend     = -25e3; % above is the fully locked depth (down dip of locked patch) (Dt)
driven      = -34e3; % driven creep, slip below this depth (this is downdip of a transition locking depth) (Dl)

% Fully locked at shallow
% lockstart   = 0e3; % below is the fully locked depth (up dip of locked patch), above creep (Dc)
% lockend     = -34e3; % above is the fully locked depth (down dip of locked patch) (Dt)
% driven      = -34e3; % driven creep, slip below this depth (this is downdip of a transition locking depth) (Dl)

% Make a box slip (for inversion test)
    % lockstart   = -0e3; % below is the fully locked depth (up dip of locked patch)
    % lockend     = -15e3; % above is the fully locked depth (down dip of locked patch)
    % driven      = -15e3; % driven creep, slip below this depth (this is downdip of a transition locking depth)
    % locked2     = -20e3;

% Set up the slip rates for the backslip and slab
sliprate = zeros(size(sig0_backslip)); 
sliprate_slab = zeros(size(sig0_slab)); 

% Set the deep slip rate and compute the stress caused by it 
sliprate(vy<=driven) = Vpl(vy<=driven); % impose creep deep slip below the locking depth
sliprate_slab(vy<=driven) = Vpl(vy<=driven); % impose creep deep slip below the locking depth
Iunlocked = (vy >= lockstart  |  (vy <= lockend & vy >= driven)); % define the unlocked patches (above lockstart and below lockend and above drive creep)

% Plot the slip ratio
f2 = figure(2);clf
    subplot(2,1,1)
        hold on
        plot(vx/1e3,Iunlocked,'k','linewidth',3)
        plot(vx/1e3,(sliprate/max(sliprate)),'Color', [0.7 0.7 0.7],'LineWidth',4) % one fault model
        plot(vx/1e3,(sliprate_slab/max(sliprate)),'--r') % slab model

        title('Slip ratio [1: unlocked & 0: locked]')
        xlabel('Distance [km]')
        ylabel('Input slip ratio')
        legend('Input locking ratio',' Input slip [onefault]','Input slip [slab]')
        grid on
        xlim([0 max(vx/1e3)])
        set(gca,'fontsize',15)


%%% Calculate the stress again if you have changed the locking rates in
%%% this script (not the mesh & slip impulses in Pylith) %%%
deepcreepstressrate = -KK * sliprate; % calculated stressing rate from backslip (or locked portion only)
KKsub = KK(Iunlocked,Iunlocked); % get the stress kernel from non-locked portion
substressrate = sig0_backslip(Iunlocked) + deepcreepstressrate(Iunlocked); % calculate stress rate from non-locked portion
substressrate_slab = sig0_slab(Iunlocked) + deepcreepstressrate(Iunlocked); % calculate stress rate from non-locked portion due to 

% Use this subset of patches to balance the stress rate
subsliprate = KKsub\substressrate;
subsliprate_slab = KKsub\substressrate_slab;

% Compute the final slip rate (give slip where it's unlocked)
sliprate(Iunlocked) = subsliprate; % note: for 3D case add strike
sliprate_slab(Iunlocked) = subsliprate_slab; 


% Plot the stressing rates
f3 = figure(3);  % plot the stressing rates like in Lindsey et al., 2021
subplot(2,1,1)
[~, lastPart] = fileparts(folder);
escapedTitle = strrep(lastPart, '_', '__');  
title(escapedTitle)
ax = gca;

    hold on
    p1out = plot(vx/1e3,sliprate/Vpl(1),'Color', 'k','LineWidth',2); % original slip imposed to model (one fault)
    p11out = plot(vx/1e3,sliprate_slab/Vpl(1),'-r','LineWidth',3); % original slip imposed to model (slab)

    % ylim([0,max(Vpl)])
    xlim([0 max(vx/1e3)])
    ylabel('Slip Deficit Rate ratio')
    xlabel('Distance [km]')
    axis tight
    grid on

    ax.XMinorGrid = 'off';  % Disable minor grid on x-axis
    ax.YMinorGrid = 'on';   % Enable minor grid on y-axis

%%% Calculate differences between backslip and ESPM %%%
% Finding the max differences between the slip rate ratio

sr_slab_ratio = sliprate_slab/Vpl(1); % ratio of slip for slab 
sr_plan_ratio = sliprate/Vpl(1); % ratio of slip for one fault

diff_ratio = abs(sr_slab_ratio-sr_plan_ratio); % absolute differences in "%"
maxdiff(counter) = max(diff_ratio);
trenchdiff(counter) = diff_ratio(1);
% lsfol(counter) = vi; % for numbers only
temp = string(split(folder,'/'));
lsfol(counter) = temp(end); % for file names

    subplot(2,1,1)
    pdif = plot(vx/1e3,diff_ratio,'--','Color',[1 1 1]*0.5);

    legend([p1out,p11out,pdif,pdif], ...
    'output slip rate [onefault]', ...
    'output slip rate [slab]', ...
    sprintf('Peak: %.1f%%, Trench: %.1f%%', maxdiff(counter)*1e2, trenchdiff(counter)*1e2), ...
    'location','northwest')
    
    set(gca,'FontSize',15)

% plot the stressing rate
subplot(2,1,2); hold on
    p2 = plot(vx/1e3,sig0_backslip,'--','Color', 0*[1 1 1],'LineWidth',2); % Stressing rate at deep creep end of plate (one fault)
    p12 = plot(vx/1e3,sig0_slab,'--c','linewidth',3); % Stressing rate at deep creep end of plate (slab)
    p3 = plot(vx/1e3,deepcreepstressrate,'Color', 0.8*[1 1 1],'LineWidth',2); % stressing rate at locked + deep creep at end of plate
    p4 = plot(vx/1e3,deepcreepstressrate+sig0_backslip,'k-','linewidth',2);%,'Color', [0.7 0.7 0.7]); % units: MPa/yr (probably?) - above without the deep creep at end of plate
    p14 = plot(vx/1e3,deepcreepstressrate+sig0_slab,'-r','linewidth',2); % units: MPa/yr (probably?) - above without the deep creep at end of plate

    xlim([0 max(vx/1e3)])
    grid on
    xlabel('Distance [km]')
    ylabel('Stressing rate')
    ylim([-5e3,0])
    legend([p2,p12,p3,p4,p14], ...
        'Sig0 [onefault]', ...
        'Sig0 [slab]', ...
        'DeepCreepStressRate', ...
        'DeepCreepStressRate+Sig0 [onefault]', ...
        'DeepCreepStressRate+Sig0 [slab]', ...
        'location','southwest')

    set(gca,'FontSize',15)


% Save figures and textfile? Textfile for "Plot_compilation_modelresults.m"
savi = 1;
if savi ~= 0
    newpath = [pwd,'/',strrep(folder, './', '')];
    fprintf('Saving figure...')

    mkdir([newpath, filesep, 'Figures'])
    saveas(f3,fullfile(newpath,'/Figures/','Slipdeficitrate_Stressingrate.png'),'png');  
    saveas(f3,fullfile(newpath,'/Figures/','Slipdeficitrate_Stressingrate.svg'),'svg');  

    %%% Save the textfiles of slip rate ratios %%%
    filename = fullfile(newpath,'SlipRateRatios.txt');  % Construct the full path

    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Write the header
    fprintf(fid, 'sr_plan_ratio\tsr_slab_ratio\tdifference_sr_ratio\tDistance[m]\tDepth[m]\n');
    
    % Write the data row-by-row
    data = [sr_plan_ratio, sr_slab_ratio, diff_ratio, vx, vy];
    fprintf(fid, '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', data');  % transpose needed to print row-by-row
    
    fclose(fid);    % Close the file

    %%% Save the peak and trench percent differences %%% 
    % Saved in the directory above the individual models
    writematrix([maxdiff' trenchdiff' lsfol'],[fileparts(newpath),'/peakdifference_sr_ratio.txt']) %

    fprintf('done\n')
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PART III: Surface displacement calculation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc

%%% Read the backslip and sig0 displacement (assuming 1 meter slip) %%%
hdf5greens = [folder,'/step_slab_greensfns-boundary.h5']; % the backslip's greens function for surface displacement
hdf5greens_sig0 = [folder_sig0,'/step_slab_greensfns_sig0-boundary.h5']; % the deep slip greens function for surface displacement

% For green's function (of locked fault), each locked fault on the displacement
cellsft_disp = h5read(hdf5greens, '/viz/topology/cells');
verticesft_disp = h5read(hdf5greens, '/geometry/vertices');
c_disp_ft = h5read(hdf5greens,'/vertex_fields/displacement'); % greens function output (assumed we used 1 meter slip to produce this displacement)

% For green's function (of ESPM, not locked fault)
cellsft_sig0_disp = h5read(hdf5greens_sig0, '/viz/topology/cells');
verticesft_sig0_disp = h5read(hdf5greens_sig0, '/geometry/vertices');
c_disp_ft_sig0 = h5read(hdf5greens_sig0,'/vertex_fields/displacement'); % greens function (assumed we used 1 meter slip to produce this displacement)


%%% For backslip %%%
% Make into a matrix
[dist,distI] = sort(verticesft_disp(1,:)/1e3); % sorting the distance order
[depth,depthI] = sort(verticesft(2,:)/1e3,'descend'); % sorting the depth order 


% Final sorting uses both sort indices
c_disp_ft_sort = c_disp_ft(:,distI,depthI); % sort the green's function due to slip impose on the locked fault
c_disp_ft_sort_sig0 = c_disp_ft_sig0(:,distI); % sig0 on locked fault in ESPM

    % Get the displacement if sliprate is from backslip
    disp_vert = (squeeze(c_disp_ft_sort(2,:,1:end-1)))*sliprate; % col is receiver, row is input slip
    disp_hori = (squeeze(c_disp_ft_sort(1,:,1:end-1)))*sliprate; % col is receiver, row is input slip
    
    % Get the displacement if sliprate is from deep creep/background (ESPM)
    disp_vert_slab = (squeeze(c_disp_ft_sort(2,:,1:end-1)))*sliprate_slab; % col is receiver, row is input slip
    disp_hori_slab = (squeeze(c_disp_ft_sort(1,:,1:end-1)))*sliprate_slab; % col is receiver, row is input slip

    % Get the displacement if there's no creep at the shallow megathrust   
    sliprate_nocreep = sliprate;
    % sliprate_nocreep(vy/1e3<=40) = 0;
    sliprate_nocreep(vx/1e3<=100) = 0;

    disp_vert_nocreep = (squeeze(c_disp_ft_sort(2,:,1:end-1)))*sliprate_nocreep; % col is receiver, row is input slip
    disp_hori_nocreep = (squeeze(c_disp_ft_sort(1,:,1:end-1)))*sliprate_nocreep; % col is receiver, row is input slip


%%% Plot the surface displacements %%%
f11 = figure(11);
    subplot(3,1,1); hold on
    % p1 = plot(dist,disp_vert,'b','linewidth',2); % vertical backslip surface displacement
    % p2 = plot(dist,disp_hori,'r','linewidth',2); % horizontal backslip surface displacement
    % 
    % p3 = plot(dist,disp_vert_slab,'c','linewidth',2); % vertical backslip surface displacement
    % p4 = plot(dist,disp_hori_slab,'m','linewidth',2); % horizontal backslip surface displacement

%%% For sig0 displacement %%%
vert_sig0 = c_disp_ft_sort_sig0(2,:)' * Vpl(1); % multiply plate rate because this is deep creep 
hori_sig0 = c_disp_ft_sort_sig0(1,:)' * Vpl(1); % multiply plate rate because this is deep creep 

    % p3 = plot(dist,vert_sig0,'--k'); % vertical sig0 surface displacement
    % p4 = plot(dist,hori_sig0,'-.k'); % horizontal sig0 surface displacement

%%% Summing both sig0 and backslip %%%
% Don't need to change the background slip because they are constant 
slab_final_vert = vert_sig0 + disp_vert; % vertical displacement final (backslip)
slab_final_hori = hori_sig0 + disp_hori; % horizontal displacement final

slab_final_vert_ESPM = vert_sig0 + disp_vert_slab; % vertical displacement final (ESPM)
slab_final_hori_ESPM = hori_sig0 + disp_hori_slab; % horizontal displacement final

slab_final_vert_nocreep = vert_sig0 + disp_vert_nocreep; % vertical displacement final (no creep)
slab_final_hori_nocreep = hori_sig0 + disp_hori_nocreep; % horizontal displacement final

    % Plot the vertical
    p11 = plot(dist,slab_final_vert,'b-','LineWidth',3); % backslip
    % p12 = plot(dist,slab_final_hori,'b-','LineWidth',3);

    p21 = plot(dist,slab_final_vert_ESPM,'r-','LineWidth',3); % ESPM
    % p22 = plot(dist,slab_final_hori_ESPM,'r-','LineWidth',3);

    % p31 = plot(dist,slab_final_vert_nocreep,'k','LineWidth',2); % no creep
    % p32 = plot(dist,slab_final_hori_nocreep,'k','LineWidth',2);

    xlim([-10 220])
    ylim([-0.015 0.02])
    set(gca,'fontsize',16)
    grid on

    ylabel({'Vertical displacement','[m/yr]'})
    xlabel('Distance from trench [km]')
    set(gca,'fontsize',15)
    title(escapedTitle)

    % plot the horizontal
    subplot(3,1,2); hold on

    % p11 = plot(dist,slab_final_vert,'b-','LineWidth',3); % backslip
    p12 = plot(dist,slab_final_hori,'b-','LineWidth',3);

    % p21 = plot(dist,slab_final_vert_ESPM,'r-','LineWidth',3); % ESPM
    p22 = plot(dist,slab_final_hori_ESPM,'r-','LineWidth',3);

    % p31 = plot(dist,slab_final_vert_nocreep,'k','LineWidth',2); % no creep
    % p32 = plot(dist,slab_final_hori_nocreep,'k','LineWidth',2);

    xlim([-10 220])
    ylim([-0.0 0.05])
    set(gca,'fontsize',16)
    grid on

    ylabel({'Horizontal displacement','[m/yr]'})
    xlabel('Distance from trench [km]')
    set(gca,'fontsize',15)
    title(escapedTitle)


% Plot the slip rate input
subplot(3,1,3); hold on
    plot(-vy/1e3,sliprate, 'b-','LineWidth',3)
    plot(-vy/1e3,sliprate_slab, 'r-','LineWidth',3)
    % plot(-vy/1e3,sliprate_nocreep, 'k-')
    
    grid on
    set(gca,'fontsize',16)

    % xlim([-10 220])
    ylabel({'Slip rate on fault','[m/yr]'})
    xlabel('Depth from trench [km]')
    title(escapedTitle)
    legend('Dislocation slip rate','ESPM slip rate','ESPM slip rate - no creep')


%%% Calculate differences between backslip and ESPM %%%
% Finding the max differences between the displacements
diff_ratio_vert = abs(slab_final_vert-slab_final_vert_ESPM); % absolute differences in their units
diff_ratio_hori = abs(slab_final_hori-slab_final_hori_ESPM); % absolute differences in their units

% Get the max %
[mdv,mdvI] = max(diff_ratio_vert); % for max
[mdh,mdhI] = max(diff_ratio_hori); % for max

% Get the trench %
tdII = find(dist == 0); % get the trench location first
[tdv,~] = max(diff_ratio_vert(tdII)); % then extract value for max difference at trench (from multiple to 1)
tdvI = find(diff_ratio_vert==tdv); % get the index of the max difference (from multiple to 1)
[tdh,~] = max(diff_ratio_hori(tdII)); % then extract value for max difference at trench
tdhI = find(diff_ratio_hori==tdh); % get the index of the max difference

% maxdiff_vert(counter) = mdv/abs(slab_final_vert(mdvI)); % %relative to BS (ESPM is % higher/lower than BS)
% maxdiff_hori(counter) = mdh/abs(slab_final_hori(mdhI)); % %relative to BS
% trenchdiff_vert(counter) = tdv/abs(slab_final_vert(tdvI)); % %find the highest differences (there might be two values at the trench)
% trenchdiff_hori(counter) = tdh/abs(slab_final_hori(tdhI)); % %find the highest differences (there might be two values at the trench)

maxdiff_vert(counter) = mdv; % %relative to BS (ESPM is % higher/lower than BS)
maxdiff_hori(counter) = mdh; % %relative to BS
trenchdiff_vert(counter) = tdv; % %find the highest differences (there might be two values at the trench)
trenchdiff_hori(counter) = tdh; % %find the highest differences (there might be two values at the trench)


f5 = figure(5); % make plot of displacement comparison between dislocation and ESPM
subplot(2,1,1); hold on; % vertical
    p11 = plot(dist,slab_final_vert,'b','LineWidth',3); % backslip
    p21 = plot(dist,slab_final_vert_ESPM,'r','LineWidth',3); % ESPM
    ylabel('Verical velocity on fault [m/yr]')

    yyaxis right
    p111 = plot(dist,diff_ratio_vert*1e3,'--','Color',[1 1 1]*0.5); % difference
    ylabel('Absolute difference [mm/yr]')
    set(gca,'ycolor',[1 1 1]*0.5) % change left y-axis color (blue)

    % p31 = plot(dist,slab_final_vert_nocreep,'k','LineWidth',2); % no creep

    xlim([-10 220])
    set(gca,'fontsize',16)
    grid on

    legend([p11,p21,p111],'Dislocation vertical disp.','ESPM vertical disp.', ...
        sprintf('Peak: %.2f mm/yr, Trench: %.2f mm/yr', maxdiff_vert(counter)*1e3, trenchdiff_vert(counter)*1e3), ...
    'location','southeast')

    title(escapedTitle)


subplot(2,1,2); hold on; % horizontal
    p12 = plot(dist,slab_final_hori,'b','LineWidth',3); % backslip
    p22 = plot(dist,slab_final_hori_ESPM,'r','LineWidth',3); % ESPM
    ylabel('Horizontal velocity on fault [m/yr]')

    yyaxis right
    p112 = plot(dist,diff_ratio_hori*1e3,'--','Color',[1 1 1]*0.5); % difference
    ylabel('Absolute difference [mm/yr]')
    set(gca,'ycolor',[1 1 1]*0.5) % change left y-axis color (blue)

    % p32 = plot(dist,slab_final_hori_nocreep,'k','LineWidth',2);

    xlim([-10 220])
    set(gca,'fontsize',16)
    grid on

    legend([p12,p22,p112],'Dislocation horizontal disp.','ESPM horizontal disp.', ...
        sprintf('Peak: %.2f mm/yr, Trench: %.2f mm/yr', maxdiff_hori(counter)*1e3, trenchdiff_hori(counter)*1e3), ...
    'location','northeast')

    xlabel('Distance from trench [km]')
    title(escapedTitle)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PART III [OPTIONAL]: Plot textfile of GNSS [if exist]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Read in Nankai's 2016 paper - note that we need to decompose their
% % magnitude and azimuth into the respective horizontals
% GPSfile = readtable('./Textfiles/NankaiTrench_GNSS_Yokota2016_reorganized_inbox.txt');
% GPS_lat = table2array(GPSfile(:,1));
% GPS_lon = table2array(GPSfile(:,2));
% GPS_EW = table2array(GPSfile(:,3));
% GPS_NS = table2array(GPSfile(:,4)); 
% GPS_tothori = table2array(GPSfile(:,5)); % [if exist - the magnitude of vector]
% GPS_angle = table2array(GPSfile(:,6)); % [if exists - the angle of vector]
% GPS_dist = table2array(GPSfile(:,7)); % Distance from trench - you need to use "get_profile.m"
% 
% % hold on to the surface displacement from the models
% figure(11)
% subplot(2,1,1); hold on
%     p15 = plot(GPS_dist,GPS_tothori/1e3,'k.','markersize',10);
% 
% pfit_tothori = polyfit(GPS_dist,GPS_tothori,3); % gnss fitting
% gps_fit_tothori = polyval(pfit_tothori,GPS_dist); hold on % total horizontal rates of GNSS
% 
%   %  p22 = plot(GPS_dist(Idist),gps_fit_tothori(Idist)/1e3,'k-'); % plot the GNSS fitting
% 
%     legend([p15],'Nankai Horizontal')


%%% Saving the figures and textfiles of surface displacement to each folder %%%
% Save figures?
savi = 1;
if savi ~= 0
    newpath = [pwd,'/',strrep(folder, './', '')];
    fprintf('Saving figure...')

    mkdir([newpath, filesep, 'Figures'])
    saveas(f5,fullfile(newpath,'/Figures/','SurfacedeformationPctDiff.png'),'png');  
    saveas(f5,fullfile(newpath,'/Figures/','SurfacedeformationPctDiff.svg'),'svg'); 

    saveas(f11,fullfile(newpath,'/Figures/','Surfacedeformation.png'),'png');  
    saveas(f11,fullfile(newpath,'/Figures/','Surfacedeformation.svg'),'svg');  

    % Save the textfiles too
    filename = fullfile(newpath,'FinalDisplacementsRatios.txt');  % Construct the full path

    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Write the header
    fprintf(fid, 'Backslip_vert[m/yr]\tBackslip_hori[m/yr]\tESPM_vert[m/yr]\tESPM_hori[m/yr]\tVertPctDiff\tHoriPctDiff\tDistance[m]\n');
    
    % Write the data row-by-row
    data = [slab_final_vert, slab_final_hori, slab_final_vert_ESPM, slab_final_hori_ESPM, diff_ratio_vert, diff_ratio_hori, dist'];
    fprintf(fid, '%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', data');  % transpose needed to print row-by-row
    
    fclose(fid);    % Close the file
    
    %%% Save the highest trench surface deformation differences between ESPM and Backslip %%% 
    % Saved in the directory above the individual models
    writematrix([maxdiff_vert' trenchdiff_vert' maxdiff_hori' trenchdiff_hori' lsfol'],[fileparts(newpath),'/peakdifference_vert_hori_meters.txt']) %
    writematrix(sliprate_slab,[(newpath),'/sliprate_ESPM.txt']) % ESPM slip rate
    writematrix(sliprate,[(newpath),'/sliprate_backslip.txt']) % Backslip sliprate (0 at locked areas)
    writematrix(sig0_backslip,[(newpath),'/sig0_backslip.txt']) % backslip sig0
    writematrix(sig0_slab,[(newpath),'/sig0_ESPM.txt']) % ESPM sig0

    fprintf('done\n')
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% % [OPTIONAL]: Plot the domain and results %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hdf5 = [folder,'/step_slab_greensfns-domain.h5'];
% cellsd = h5read(hdf5, '/viz/topology/cells');
% verticesd = h5read(hdf5, '/geometry/vertices');
% 
% plot_domain_pylith(cellsd,verticesd,'r')
% plot_domain_pylith(cellsgs,verticesgs)
% plot_domain_pylith_quiver(cellsgs,verticesgs,vertexftgs)
% 
% xsver = verticesft(1,:);
% ysver = verticesft(2,:);

% Plot the displacement (slip)
% for ii = 1:size(c_slip_ft_sort,3)
% 
%     sliporder(ii,:) = c_slip_ft_sort(2,:,ii); % see the order of the patch slip imposed
% 
%     figure(1),hold on
% 
%         quiver(sortx/1e3,sorty/1e3,c_slip_ft_sort(1,:,ii),c_slip_ft_sort(2,:,ii)) % 2nd slip for thrust, 1st for strike-slip
% 
% end
% 
%     xlim([-50 500])
%     ylim([-100 10])
%     xlabel('Distance (km)')
%     ylabel('Depth (km)')
%     grid on
    % axis square


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% % [OPTIONAL]: Plot the strain/stress (traction) at every slip 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f3 = figure(3); clf,hold on
% 
% for ii = 1:size(c_trac_ft_sort,3)
% 
%     % Extract the components (Traction) [component,location,time]
%     tno = c_trac_ft_sort(1,:,ii); % normal 
%     tas = c_trac_ft_sort(2,:,ii); % shear
% 
%     tno_sort = c_trac_ft_sort(1,locked_I_sortreceiver,ii); 
%     tas_sort = c_trac_ft_sort(2,locked_I_sortreceiver,ii); 
% 
%         subplot(2,1,1)
%         % plot(sortx(1:108)/1e3,tas_sort(1:108)); hold on % make sure you get the locked portion only or your choice
%         plot(sortx_locked/1e3,tas_sort); hold on % make sure you get the locked portion only or your choice
%         ylabel('Shear')
%         xlabel('Distance')
%         set(gca,'fontsize',20)
% 
%         grid on
% 
%         subplot(2,1,2)
%         % plot(sortx(1:108)/1e3,tno_sort(1:108)); hold on
%         plot(sortx_locked/1e3,tno_sort); hold on
%         ylabel('Normal')
%         xlabel('Distance')
% 
%         % xlim([0 500])
%         % ylim([0 40])
% end
% 
% title('Sorted by distance from sorted x-axis')
% 
%     grid on
%     set(gca,'fontsize',20)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Plot the matrix of slip and GNSS %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f4 = figure(4);hold on
%     % pcolor(squeeze(vertexftgs(1,locked_I_sortreceiver,:))),shading flat
%     imagesc(squeeze(c_slip_ft_sort(2,:,:)))
% 
%     axis tight
%     colorbar
% 
%     xlabel('Slip on patch number - Check sorting here')
%     ylabel('Observers number (GNSS/fault patch)')
%     set(gca,'fontsize',20)
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Plot the stress kernel matrix (like green's function)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% f5 = figure(5);hold on
% f5.Position = [1           1        1920        1096];
% subplot(2,2,1)
%     pcolor(squeeze(c_trac_ft_sort(2,:,:))),shading flat
% 
%     axis square
%     colorbar
% 
%     xlabel('Patch number (source)')
%     ylabel('Patch number (reciever)')
%     title('Shear stress interaction kernel, sorted by distance x1')
% 
%     set(gca,'fontsize',20)
% 
% subplot(2,2,2)
%     pcolor(squeeze(c_trac_ft_sort(1,:,:))),shading flat
% 
%     axis square
%     colorbar
% 
%     xlabel('Patch number (source)')
%     ylabel('Patch number (reciever)')
%     title('Normal stress interaction kernel, sorted by distance x1')
% 
%     set(gca,'fontsize',20)
% 
% subplot(2,2,3)
%     pcolor(squeeze(log(abs(c_trac_ft_sort(2,:,:))))+1),shading flat
% 
%     axis square
%     colorbar
% 
%     xlabel('Patch number (source)')
%     ylabel('Patch number (reciever)')
%     title('Shear stress interaction kernel, log(abs(traction))+1')
% 
%     set(gca,'fontsize',20)
% 
% subplot(2,2,4)
%     pcolor(squeeze(log(abs(c_trac_ft_sort(1,:,:))))+1),shading flat
% 
%     axis square
%     colorbar
% 
%     xlabel('Patch number (source)')
%     ylabel('Patch number (reciever)')
%     title('Normal stress interaction kernel, log(abs(traction))+1')
% 
%     set(gca,'fontsize',20)
% 
% 
% % Save figure?
% savi = 0;
% if savi == 1
%     newpath = [pwd,'/',strrep(folder, './', '')];
%     fprintf('Saving figure...')
%     mkdir([newpath, filesep, 'Figures'])
%     saveas(f3,fullfile(newpath,'/Figures/','StressPatches_sorted.png'),'png');  
%     saveas(f4,fullfile(newpath,'/Figures/','SlipGNSS_matrix_sorted.png'),'png');  
%     saveas(f5,fullfile(newpath,'/Figures/','StressKernels_sorted.png'),'png');  
% 
%     fprintf('Saving stress kernel textfiles...')
%     mkdir([newpath, filesep, 'Textfiles'])
%     filename = fullfile(newpath, 'Textfiles', 'ShearStressKernel.txt');  % Construct the full path
%     writematrix(c_trac_ft_sort(1,:,:), filename);
%     filename = fullfile(newpath, 'Textfiles', 'NormalStressKernel.txt');  % Construct the full path
%     writematrix(c_trac_ft_sort(2,:,:), filename);
%     fprintf('done\n')
% end