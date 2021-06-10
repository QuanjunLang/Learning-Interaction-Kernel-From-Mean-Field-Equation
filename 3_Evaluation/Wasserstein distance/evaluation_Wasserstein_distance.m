function ws_dist = evaluation_Wasserstein_distance(sysInfo, infer_output, U, new_initial)
% compute the Wasserstein distance between the data U and the data
% gernerated by the estimated interaction kernel

%% Use estimated potential to generate new solution u
% fprintf('Try different initial values ... ')
% new_initial = 'DN_2_1';
% diff_ini_plotON = 0;
% get_u_from_new_gradK;
% fprintf('Done!\n')

% Use estimated potential to genreate solution u


%% Load parameters
xgrid=  sysInfo.xgrid;
dt = sysInfo.dt;
T = sysInfo.T;


phi_kernel = infer_output.phi_kernel;
%% Make a function handle of the estimated interaction kernel
temp = fit(xgrid', phi_kernel.est ,'smoothingspline');
phi_kernel_function_handle.est = @(x) reshape(temp(x), size(x));
phi_kernel_function_handle.true = sysInfo.phi_kernel;

%% old initial, true kernel
U1 = U;  

%% old initial, estimated kernel
sysInfo.phi_kernel = phi_kernel_function_handle.est;
U2 = SPCC_Imp_solver(sysInfo, 0);

%% new initial, true kernel
sysInfo.phi_kernel = phi_kernel_function_handle.true;
sysInfo.Initial = new_initial;
sysInfo = get_sysInfo_details(sysInfo);

U3 = SPCC_Imp_solver(sysInfo, 0);

%% new initial, estimated kernel
sysInfo.Initial = new_initial;
sysInfo = get_sysInfo_details(sysInfo);
sysInfo.phi_kernel = phi_kernel_function_handle.est;

U4 = SPCC_Imp_solver(sysInfo, 0);


%% compute distance
ws_dist.old_initial = zeros(sysInfo.TN+1, 1);
ws_dist.new_initial = zeros(sysInfo.TN+1, 1);

for i = 1:sysInfo.TN
    ws_dist.old_initial(i) = ws_distance(U1(:, i), U2(:, i), 2);
    ws_dist.new_initial(i) = ws_distance(U3(:, i), U4(:, i), 2);
end









end