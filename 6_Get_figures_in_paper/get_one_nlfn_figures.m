% Givern the density u(x,t)
% Infer the interaction kernel phi

%%
% clear all;
% close all;
% clc
% all_dir = add_path_create_dir();


display_sysInfo(sysInfo)
 
%% load specified parameters
% sysInfo = settings_DEs();
% infer = settings_inference();

%% sysInfo settings (optional)
% sysInfo.L = 10;
% sysInfo.M = 300;
% sysInfo.T = 1;
% sysInfo.dt = 0.001;

% sysInfo.Initial = 'DiffN_2_0.5_-3_1';
% sysInfo.nlfn = 'AB_2_-0.5_cut_threshold_0.05';
% sysInfo.v = 0.01;

% sysInfo.Initial = 'DN_1_0.5';
% sysInfo.nlfn = 'power_3';
% sysInfo.v = 1;

% sysInfo.v = 0.1; 
% sysInfo.nlfn = 'Opinion_P1_3_-1_P2_4_2';
% sysInfo.Initial = 'DiffN_-2_1_-4_0.5_2_1';
% 
% sysInfo = get_sysInfo_details(sysInfo);



% settings_AB;

%% Generate/Load Data U
fprintf(['Generate data U with M = ', num2str(sysInfo.M),'...\n'])
U = generate_data(sysInfo, all_dir, 1);
plot_U(U, sysInfo, all_dir, 1);




%% Inference Computational kernel
fprintf('Generate integral Kernel G...\n')
G = generate_integration_kernel(sysInfo, U, all_dir);


%% Inference settings (optional)
% infer.type = 'spline';
% infer.basis_num = 'deg_2_knotnum_26';
% infer.b_method = 'ptgrad';
% infer.c_method = 'Tik_H1_unif_optimal_pinv';
% infer.rho_bdry = 5e-4;
% infer = get_infer_details(infer);

%% Inference with spline (optimal dimension)
% spline_candidate_N = 30;
[infer_output, infer] = inference_kernel_optimal_dimension(U, sysInfo, infer, G, spline_candidate_N, all_dir, 0);
fprintf('optimal dimension for spline basis is %d\n', infer.n);
%% Inference with RKHS (optimal dimension)
% rkhs_infer = infer;
% rkhs_infer.type = 'rkhs_rho_mesh+deg_2_knotnum_100';
% rkhs_infer.basis_num = 'n_40';
% rkhs_infer.c_method = 'Tik_I_unif_optimal_pinv';
% rkhs_infer = get_infer_details(rkhs_infer);
% 
% rkhs_candidate_N = 40;
[infer_output_rkhs, rkhs_infer] = inference_kernel_optimal_dimension(U, sysInfo, rkhs_infer, G, rkhs_candidate_N, all_dir, 0);
fprintf('optimal dimension for rkhs basis is %d\n', rkhs_infer.n);


plot_2_inference(sysInfo, infer, infer_output, infer_output_rkhs, G, all_dir, 1)
%% Wasserstein distance
new_initial = 'DN_2_1';
ws_dist = evaluation_Wasserstein_distance(sysInfo, infer_output, U, new_initial);
plot_Wasserstein_distance(ws_dist, sysInfo, all_dir, 1);

%% Free energy
free_Engy = evaluation_Free_energy(sysInfo, infer_output, U);
plot_Free_energy(free_Engy, sysInfo, all_dir, 1)

%% meshsize compare
all_M = [30, 40, 50, 60, 100, 120, 150, 200, 250, 300];
big_M = 3000;
[Full_U, all_infer_output, all_Error, all_G] = meshsize_compare(sysInfo, infer, all_dir, all_M, big_M);
plot_meshsize_compare(sysInfo, all_infer_output, all_M, all_dir, 1);