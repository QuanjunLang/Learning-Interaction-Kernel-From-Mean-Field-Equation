% Generate all figures in the paper
%%
clear all;
close all;
clc
all_dir = add_path_create_dir();

%% L curve and dimension selection
L_D_plot(all_dir)
%% load specified parameters
sysInfo = settings_DEs(1);
infer = settings_inference();

%% cubic potential
fprintf('cubic potential...\n')
settings_power_3;
get_one_nlfn_figures;

%% Opinion dynamics
fprintf('\nOpinion dynamics potential...\n')
settings_Opinion;
get_one_nlfn_figures;

%% AB potential
fprintf('\nRepulsive aggresive potential...\n')
settings_AB;
get_one_nlfn_figures;
