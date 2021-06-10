%% sysInfo settings (optional)
sysInfo.L = 10;
sysInfo.M = 300;
sysInfo.T = 1;
sysInfo.dt = 0.001;

sysInfo.Initial = 'DN_1_0.5';
sysInfo.nlfn = 'power_3';
sysInfo.v = 1;

sysInfo = get_sysInfo_details(sysInfo);

%% inference settings (spline)
infer.type = 'spline';
infer.basis_num = 'deg_2_knotnum_26';
infer.b_method = 'ptgrad';
infer.c_method = 'Tik_H1_unif_optimal_pinv';
infer.rho_bdry = 5e-4;
infer = get_infer_details(infer);

spline_candidate_N = 30;
%% inference settings (rkhs)
rkhs_infer = infer;
rkhs_infer.type = 'rkhs_rho_mesh+deg_2_knotnum_100';
rkhs_infer.basis_num = 'n_40';
rkhs_infer.c_method = 'Tik_I_unif_optimal_pinv';
rkhs_infer = get_infer_details(rkhs_infer);

rkhs_candidate_N = 40;
