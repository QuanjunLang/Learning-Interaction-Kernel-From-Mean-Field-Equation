%% sysInfo settings (optional)
sysInfo.L = 10;
sysInfo.M = 300;
sysInfo.T = 1;
sysInfo.dt = 0.001;

sysInfo.Initial = 'DiffN_2_0.5_-3_1';
sysInfo.nlfn = 'AB_2_-0.5_cut_threshold_0.05';
sysInfo.v = 0.01;

sysInfo = get_sysInfo_details(sysInfo);
%% inference settings (spline)
infer.type = 'spline';
infer.basis_num = 'deg_1_knotnum_10';
infer.b_method = 'ptgrad';
infer.c_method = 'Tik_H1_unif_optimal_pinv';
infer.rho_bdry = 5e-4;
infer = get_infer_details(infer);

spline_candidate_N = 30;
%% inference settings (rkhs)
rkhs_infer = infer;
rkhs_infer.type = 'rkhs_rho_mesh+deg_1_knotnum_100';
rkhs_infer.basis_num = 'n_10';
rkhs_infer.c_method = 'Tik_I_unif_optimal_pinv';
rkhs_infer = get_infer_details(rkhs_infer);

rkhs_candidate_N = 40;
