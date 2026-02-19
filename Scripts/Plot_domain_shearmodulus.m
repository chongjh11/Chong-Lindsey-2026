%% Plot domain shear modulus 
%
%   NOTES:
%   - This plots the shear modulus/bulk modulus produced from Pylith
%
%   INSTRUCTIONS:
%   1. Ensure the right path to the folder
%   2. Check "domainfiles" if you have the correct domains with material
%   properties, change accordingly
%
%   Last modified on 17-Sep-2025
%   by Jeng Hann, Chong
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%%% List the info files - CHANGE ACCORDINGLY %%%
% You might need to remove the "accretionary wedge" if your mesh doesn't have it
% domainfiles = ["/step_slab_greensfns-continental_mantle_info.h5"; "/step_slab_greensfns-oceanic_mantle_info.h5";"/step_slab_greensfns-ocean_crust_info.h5";"step_slab_greensfns-accretionary_wedge_info.h5"];
domainfiles = ["/step_slab_greensfns-continental_mantle_info.h5"; "/step_slab_greensfns-oceanic_mantle_info.h5";"/step_slab_greensfns-ocean_crust_info.h5"];

shrms_comp2 = []; % leave empty here
cellscm = []; % leave empty

allfolds = readtable('./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/30deg60thc/listmatwed30deg.txt', 'ReadVariableNames', false,'Delimiter','.');
% allfolds = readtable('./Deformation/Models_final_slab_thickness_2d/listmesh.txt', 'ReadVariableNames', false,'Delimiter','.');


for vi = 1:size(allfolds,1) % Change accordingly - this to get the folder versions

close all
clc

    for di = 1:length(domainfiles) % loop for plotting each domain
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART I: COMPILE AND EXTRACT DATA FOR COLOR 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Load Input %%% CHANGE ACCORDINGLY
    % folder_sig0 = ['./Deformation/Models_final_slab_wedge_elastic_2d/1GPaWedge_30GPaSlab_30GPaCrust_70GPaMantle/elastic_slab_sig0_v',num2str(vi)]; 
    % folder = ['./Deformation/Models_final_slab_wedge_elastic_2d/1GPaWedge_30GPaSlab_30GPaCrust_70GPaMantle/elastic_slab_v',num2str(vi)];
    % folder_sig0 = ['./Deformation/Models_final_slab_wedge_elastic_2d/test/elastic_slab_sig0_v',num2str(vi)]; 
    % folder = ['./Deformation/Models_final_slab_wedge_elastic_2d/test/elastic_slab_v',num2str(vi)];
    % folder_sig0 = ['./Deformation/Models_final_slab_elastic_2d/test/elastic_slab_sig0_v',num2str(vi)]; 
    % folder = ['./Deformation/Models_final_slab_elastic_2d/test/elastic_slab_v',num2str(vi)] ;
    % folder_sig0 = ['./Deformation/Models_final_slab_matwedge_elastic_2d/10deg40thc/elastic_slab_sig0_v',num2str(vi)]; 
    % folder = ['./Deformation/Models_final_slab_matwedge_elastic_2d/10deg40thc/elastic_slab_v',num2str(vi)];

    naming = char(table2cell(allfolds(vi,1))); % needed for naming the folders
    folder_sig0 = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/30deg60thc/elastic_slab_sig0_',naming];
    folder = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/Forward/30deg60thc/elastic_slab_',naming];
    % folder_sig0 = ['./Deformation/Models_final_slab_thickness_2d/elastic_slab_sig0_',naming];
    % folder = ['./Deformation/Models_final_slab_thickness_2d/elastic_slab_',naming];

    % folder_sig0 = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/elastic_slab_sig0']
    % folder = ['./Deformation/Models_final_slab_varymatwedge_elastic_2d/elastic_slab']

    % Extract the properties (don't need to change)
    hdf5cm = append(folder,'/',domainfiles(di));
    cellscm = h5read(hdf5cm, '/viz/topology/cells');
    verticescm = h5read(hdf5cm, '/geometry/vertices');
    shrmodcm = h5read(hdf5cm,'/vertex_fields/shear_modulus');
    blkmodcm = h5read(hdf5cm,'/vertex_fields/bulk_modulus');
    
     % plot_domain_pylith(cellscm,verticescm,'r') % plot the entire domain
    
    % round to nearest GPa (for making the colors smoother)
    shrmodcm = round(shrmodcm/1e9);

    % Compiled information (might not use everything)
    cells_comp{di,1} = cellscm;
    verts_comp{di,1} = verticescm;
    shrms_comp{di,1} = shrmodcm;
    shrms_comp2(length(shrms_comp2)+1:length(shrmodcm)+length(shrms_comp2)) = shrmodcm;
    blkms_comp{di,1} = blkmodcm;

    end

    %%% Make colors based on a scale of your values (E.g., GPa) %%%
    % unishrmod = unique(shrms_comp2); % based on whatever your model is
    unishrmod = 1:70; % give fixed values
    cmap = turbo(length(unishrmod)); % recheck this to see if the plot really took the values correctly
    idx_comp = [];

    % Match the domains with the colors in the cell format
    for di = 1:length(domainfiles)   % since you fixed it

        [~, idx] = ismember(shrms_comp{di,1}, unishrmod);
        cls_comp{di,1} = cmap(idx, :); % RGB
        idx_comp = [idx_comp;idx']; % index for RGB to shear mod in whole matrix
        idx_comp_2{di,1} = idx'; % cells version 

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 2: RELOAD DOMAINS AND PLOT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    shrmodnow = [];
    vertices_x_now = [];
    vertices_y_now = [];
    cells_x_now = [];
    cells_y_now = [];
    
    for di = 1:length(domainfiles)

    shrmodnow(length(shrmodnow)+1:length(shrmodnow)+length(shrms_comp{di,1}),1) = shrms_comp{di,1}'; % reclassify the domains
    vertices_x_now(length(vertices_x_now)+1:length(vertices_x_now)+length(verts_comp{di,1}),1) = verts_comp{di,1}(1,:)'; % xaxis - reclassify vertices for each domain
    vertices_y_now(length(vertices_y_now)+1:length(vertices_y_now)+length(verts_comp{di,1}),1) = verts_comp{di,1}(2,:)'; % depth - reclassify vertices for each domain
    cells_x_now(length(cells_x_now)+1:length(cells_x_now)+length(cells_comp{di,1}),1) = cells_comp{di,1}(1,:)'; % xaxis - reclassify vertices for each domain
    cells_y_now(length(cells_y_now)+1:length(cells_y_now)+length(cells_comp{di,1}),1) = cells_comp{di,1}(2,:)'; % depth - reclassify vertices for each domain

    end

    %%% Do combined vertices plot %%%
    % Preallocate empty arrays
    clear V_all F_all C_all S_all

    V_all = [];
    F_all = [];
    C_all = [];
    S_all = [];
    pdi = [];

    for di = 1:length(domainfiles)

        V = verts_comp{di,1}.';
        F = cells_comp{di,1}.';
        C = cls_comp{di,1};
        S = shrms_comp{di,1}.';

        % Example: if F has negative values, convert to positive indices
        % Map negative indices to positive
        minF = min(F(:));
        if minF <= 0
            F_corrected = F - minF + 1;  % shift all indices so the smallest is 1
        end
 
        % Offset F by the number of vertices already in V_all
        V_all = [V_all; V];
        F_all = [F_all; F_corrected];
        C_all = [C_all; C];
        S_all = [S_all; S];


        % Reorganize the colors (based on all the vertices)
        idx2 = idx_comp_2{di,1};
        uniqueColors = cls_comp{di,1};

        f2 = figure(2);
        f2.Position = [297 888 1039 371];

        patch('Faces', F_corrected, ...
            'Vertices', V, ...
            'FaceVertexCData', idx2, ...   % RGB per face
            'FaceColor', 'flat', ...    % 'flat' for face-wise coloring
            'CDataMapping', 'direct', ... % use direct mapping (ignore colormap scaling)
            'EdgeColor', 'none', ...
            'LineWidth', 0.2)

        hold on
        colorbar

    end
        colormap(cmap)

        xlim([-1 1]*1e5)
        ylim([-100 10]*1e3)
    
        grid on
        [~, lastPart] = fileparts(folder);
        escapedTitle = strrep(lastPart, '_', '__');  
        title(escapedTitle)
        xlabel('Distance [km]');ylabel('Depth [km]')
        set(gca,'fontsize',15)



    %%% Do vertex only plot %%%
    f1 = figure(1); hold on;
    f1.Position = [297 288 1039 371];
    colormap(turbo)

    p1 = scatter(vertices_x_now,vertices_y_now,40, shrmodnow,'filled'); % plot row1 vs row2 for each column

    colorbar
    clim([0 max(shrmodnow)]) % start from zero because other models have 1e9 rather than 1e10

    xlim([-2 2]*1e5)
    ylim([-100 10]*1e3)

    grid on
    [~, lastPart] = fileparts(folder);
    escapedTitle = strrep(lastPart, '_', '__');  
    title(escapedTitle)
    xlabel('Distance [km]');ylabel('Depth [km]')
    set(gca,'fontsize',15)

    fprintf('done plotting\n')

    % OPTIONAL for saving
    savi = 1;
    if savi ~= 0
        newpath = [pwd,'/',strrep(folder, './', '')];
        fprintf('Saving figure...')
    
        mkdir([newpath, filesep, 'Figures'])
        saveas(f2,fullfile(newpath,'/Figures/','ShearModulusDomain.png'),'png');  
        saveas(f2,fullfile(newpath,'/Figures/','ShearModulusDomain.svg'),'svg'); 
    end
end