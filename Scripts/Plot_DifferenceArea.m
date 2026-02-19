%% Plot the difference in area between the inverse and forward models
% 
%   - Get the compiled text files list of area differences textfile from "ForwardInverse_2d_loop.m"
%   - Change the file name accordingly, it will produce the area and % area
%   from the textfile you gave
%
%   Last modified on 27-Nov-2025
%   by Jeng Hann Chong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

clear all
close all
clc

% OPTIONAL - only plot select models (choose accordingly)
slc = [1:5:50];

% Load the input
a = readtable('./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/20deg40thc/area_difference_inv_fwd_fulloffonshore.txt','Delimiter',',');
a = a(slc,:);

mean_area_diff_fulloffonshore = table2array(a(:,1));
max_area_diff_fulloffonshore = table2array(a(:,2));
min_area_diff_fulloffonshore = table2array(a(:,3));
mean_pct_diff_fulloffonshore = table2array(a(:,4));
max_pct_diff_fulloffonshore = table2array(a(:,5));
min_pct_diff_fulloffonshore = table2array(a(:,6));
names = a(:,end);
filenum = (1:size(a,1))';

b = readtable('./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/20deg40thc/area_difference_inv_fwd_fullonshore.txt','Delimiter',',');
b = b(slc,:);

mean_area_diff_fullonshore = table2array(b(:,1));
max_area_diff_fullonshore = table2array(b(:,2));
min_area_diff_fullonshore = table2array(b(:,3));
mean_pct_diff_fullonshore = table2array(b(:,4));
max_pct_diff_fullonshore = table2array(b(:,5));
min_pct_diff_fullonshore = table2array(b(:,6));

c = readtable('./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/20deg40thc/area_difference_inv_fwd_horioffonshore.txt','Delimiter',',');
c = c(slc,:);

mean_area_diff_horioffonshore = table2array(c(:,1));
max_area_diff_horioffonshore = table2array(c(:,2));
min_area_diff_horioffonshore = table2array(c(:,3));
mean_pct_diff_horioffonshore = table2array(c(:,4));
max_pct_diff_horioffonshore = table2array(c(:,5));
min_pct_diff_horioffonshore = table2array(c(:,6));

d = readtable('./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/20deg40thc/area_difference_inv_fwd_horionshore.txt','Delimiter',',');
d = d(slc,:);

mean_area_diff_horionshore = table2array(d(:,1));
max_area_diff_horionshore = table2array(d(:,2));
min_area_diff_horionshore = table2array(d(:,3));
mean_pct_diff_horionshore = table2array(d(:,4));
max_pct_diff_horionshore = table2array(d(:,5));
min_pct_diff_horionshore = table2array(d(:,6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot absolute area difference %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure(1);
s1 = subplot(3,1,1); hold on
    
    fill([filenum; flipud(filenum)],[min_area_diff_fulloffonshore; flipud(max_area_diff_fulloffonshore)],[0.8 0.1 0.1], 'EdgeColor','none','FaceAlpha',0.5) 
    p1 = plot(mean_area_diff_fulloffonshore,'r.-');

    fill([filenum; flipud(filenum)],[min_area_diff_horioffonshore; flipud(max_area_diff_horioffonshore)],[0.1 0.1 0.8], 'EdgeColor','none','FaceAlpha',0.5) 
    p2 = plot(mean_area_diff_horioffonshore,'b.-');

    title('Off & onshore data')
    legend([p1,p2],'Vert & Hori','Hori','location','southeast')
    grid on

s2 = subplot(3,1,2); hold on

    fill([filenum; flipud(filenum)],[min_area_diff_fullonshore; flipud(max_area_diff_fullonshore)],[0.8 0.1 0.1], 'EdgeColor','none','FaceAlpha',0.5) 
    p3 = plot(mean_area_diff_fullonshore,'r.-');

    fill([filenum; flipud(filenum)],[min_area_diff_horionshore; flipud(max_area_diff_horionshore)],[0.1 0.1 0.8], 'EdgeColor','none','FaceAlpha',0.5) 
    p4 = plot(mean_area_diff_horionshore,'b.-');

    title('Onshore data only')
    legend([p3,p4],'Vert & Hori','Hori','location','southeast')
    grid on

    han=axes(f1,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';

    ylabel('Area difference between input and predicted')
    xlabel('Model number')
    fontsize(gcf,14,'points')
    
    set(gca,'fontsize',16)
    sgtitle('Negative: predict > input; Positive: predict < input')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot percent (%) area difference %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = figure(2);
s21 = subplot(3,1,1); hold on
    
    fill([filenum; flipud(filenum)],[min_pct_diff_fulloffonshore; flipud(max_pct_diff_fulloffonshore)]*1e2,[0.8 0.1 0.1], 'EdgeColor','none','FaceAlpha',0.5) 
    p21 = plot(mean_pct_diff_fulloffonshore*1e2,'r.-');

    fill([filenum; flipud(filenum)],[min_pct_diff_horioffonshore; flipud(max_pct_diff_horioffonshore)]*1e2,[0.1 0.1 0.8], 'EdgeColor','none','FaceAlpha',0.5) 
    p22 = plot(mean_pct_diff_horioffonshore*1e2,'b.-');

    title('Off & onshore data')
    legend([p21,p22],'Vert & Hori','Hori','location','southeast')
    grid on

    ylim([-50 30]) % set the ylim here
    ax = gca;
    ax.YMinorGrid = 'on';

s22 = subplot(3,1,2); hold on

    fill([filenum; flipud(filenum)],[min_pct_diff_fullonshore; flipud(max_pct_diff_fullonshore)]*1e2,[0.8 0.1 0.1], 'EdgeColor','none','FaceAlpha',0.5) 
    p23 = plot(mean_pct_diff_fullonshore*1e2,'r.-');

    fill([filenum; flipud(filenum)],[min_pct_diff_horionshore; flipud(max_pct_diff_horionshore)]*1e2,[0.1 0.1 0.8], 'EdgeColor','none','FaceAlpha',0.5) 
    p24 = plot(mean_pct_diff_horionshore*1e2,'b.-');

    ylim([-50 30]) % set the ylim here
    ax = gca;
    ax.YMinorGrid = 'on';

    title('Onshore data only')
    legend([p23,p24],'Vert & Hori','Hori','location','southeast')
    grid on

    han=axes(f2,'visible','off'); 
    han.Title.Visible='on';
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';

    ylabel('Area difference between input and predicted [%]')
    xlabel('Model number')
    fontsize(gcf,14,'points')
    
    set(gca,'fontsize',16)
    sgtitle('Negative: predict > input; Positive: predict < input')

