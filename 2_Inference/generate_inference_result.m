function infer_output = generate_inference_result(infer, sysInfo, U, G, plot_L_curve, all_dir, big_M)

%% get file name
%% 
if nargin == 6
    infer_output_name = [all_dir.infer_dir, 'nlfn_',sysInfo.nlfn, '_Initial_',sysInfo.Initial];
    infer_output_name = [infer_output_name, '_M_', num2str(sysInfo.M), '_L_', num2str(sysInfo.L)];
    infer_output_name = [infer_output_name, '_dt_',num2str(sysInfo.dt),'_v_',num2str(sysInfo.v),'_',infer.type, '_', infer.basis_num, '_', infer.c_method, '_', infer.b_method, '.mat'];
elseif nargin == 7
    infer_output_name = [all_dir.infer_dir, 'nlfn_',sysInfo.nlfn, '_Initial_',sysInfo.Initial];
    infer_output_name = [infer_output_name, '_bigM_', num2str(big_M), '_sparseM_', num2str(sysInfo.M), '_L_', num2str(sysInfo.L)];
    infer_output_name = [infer_output_name, '_dt_',num2str(sysInfo.dt),'_v_',num2str(sysInfo.v),'_',infer.type, '_', infer.basis_num, '_', infer.c_method, '_', infer.b_method, '.mat'];
end

if ~exist(infer_output_name, 'file')
    infer_output = inference_kernel_all(infer, sysInfo, U, G, plot_L_curve); 
    save(infer_output_name, 'infer_output');
else
    load(infer_output_name, 'infer_output');
end
end