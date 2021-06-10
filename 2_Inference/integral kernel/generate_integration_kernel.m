function G = generate_integration_kernel(sysInfo, U, all_dir, big_M)
% generate kernel
% If exists, then load kernel. If not, then generate

%% 
if nargin == 3
    kernel_name = [all_dir.infer_dir, 'nlfn_',sysInfo.nlfn, '_Initial_',sysInfo.Initial];
    kernel_name = [kernel_name, '_M_', num2str(sysInfo.M), '_L_', num2str(sysInfo.L)];
    kernel_name = [kernel_name, '_dt_',num2str(sysInfo.dt),'_v_',num2str(sysInfo.v),'_kernel.mat'];
else
    kernel_name = [all_dir.infer_dir, 'nlfn_',sysInfo.nlfn, '_Initial_',sysInfo.Initial];
    kernel_name = [kernel_name, '_bigM_', num2str(big_M), '_sparseM_', num2str(sysInfo.M), '_L_', num2str(sysInfo.L)];
    kernel_name = [kernel_name, '_dt_',num2str(sysInfo.dt),'_v_',num2str(sysInfo.v),'_kernel.mat'];
end

if ~exist(kernel_name, 'file')
    G = inference_get_integration_kernel(sysInfo, U);
    save(kernel_name, 'G');
else
    load(kernel_name, 'G');
end


end