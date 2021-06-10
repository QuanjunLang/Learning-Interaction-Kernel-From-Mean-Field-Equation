% Givern the density u(x,t)
% Infer the interaction kernel phi

%%
clear all;
close all;
clc
all_dir = add_path_create_dir();
 
%% load specified parameters
sysInfo = settings_DEs(1);
infer = settings_inference();
display_sysInfo(sysInfo);
% %% sysInfo settings (optional)
% sysInfo.L = 10;
% sysInfo.M = 300;
% sysInfo.T = 1;
% sysInfo.dt = 0.001;
% sysInfo.Initial = 'DN_1_0.5';
% sysInfo.nlfn = 'power_3';
% sysInfo.v = 1;
% sysInfo = get_sysInfo_details(sysInfo);

% %% Inference settings (optional)
% infer.type = 'spline';
% infer.basis_num = 'deg_2_knotnum_26';
% infer.b_method = 'ptgrad';
% infer.c_method = 'Tik_H1_unif_optimal_pinv';
% infer.rho_bdry = 5e-4;
% infer = get_infer_details(infer);


%% Generate/Load Data U
fprintf(['Generate data U with M = ', num2str(sysInfo.M),'...\n'])
U = generate_data(sysInfo, all_dir, 1);
plot_U(U, sysInfo, all_dir, 0);

%% Inference Computational kernel
fprintf('Generate integral Kernel G...\n')
G = generate_integration_kernel(sysInfo, U, all_dir);

%% Inference (optimal dimension)

% 1. use specified dimension

% plot_L_curve = 0;
% infer_output = inference_kernel_all(infer, sysInfo, U, G, plot_L_curve);                   

% 2. find and use optimal dimension
candidate_N = 30;
[infer_output, infer] = inference_kernel_optimal_dimension(U, sysInfo, infer, G, candidate_N, all_dir, 0);
fprintf('optimal dimension is %d\n', infer.n);


Error = evaluation_error(sysInfo, infer_output, G);
plot_inference(sysInfo, infer, infer_output, G, all_dir, 0)


%% Wasserstein distance
new_initial = 'DN_2_1';
ws_dist = evaluation_Wasserstein_distance(sysInfo, infer_output, U, new_initial);
plot_Wasserstein_distance(ws_dist, sysInfo, all_dir, 0);

%% Free energy
free_Engy = evaluation_Free_energy(sysInfo, infer_output, U);
plot_Free_energy(free_Engy, sysInfo, all_dir, 0)


