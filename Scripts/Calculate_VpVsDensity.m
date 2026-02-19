%% Calculate Vp,Vs given shear mod, density and poisson ratio

clear all
close all
clc

%%% Input your desired shear modulus and others %%%
shrmod = (1:10)*1e9; % shear modulus (Pascal)
pr = 0.25; % poisson ratio [unitless]
dens = 2650; % density [kg/m3]

% Vs
Vs = sqrt(shrmod/dens); % m/s

% Vp
Vp = sqrt((Vs.^2)*(2*(1-pr))/(1-2*pr)); % m/s

% Bulk mod
% K = (2*shrmod*(1+pr))/(3*(1-2*pr)); % pascal
K = dens*(Vp.^2 - (4/3)*Vs.^2);

A = [shrmod',Vs',Vp'];

% Compile
compiled = [Vs; Vp; shrmod];
