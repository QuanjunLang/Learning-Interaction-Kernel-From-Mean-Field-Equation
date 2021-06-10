function U = generate_data(sysInfo, all_dir, progressON)
% generate data
% If exists, then load data. If not, then generate

%%
data_name = [all_dir.data_dir, 'nlfn_',sysInfo.nlfn, '_Initial_',sysInfo.Initial];
data_name = [data_name, '_M_', num2str(sysInfo.M), '_L_', num2str(sysInfo.L), '_dt_',num2str(sysInfo.dt),'_v_',num2str(sysInfo.v),'.mat'];

if ~exist(data_name, 'file')
    U = SPCC_Imp_solver(sysInfo, progressON);
    save(data_name, 'U');
else
    load(data_name, 'U');
end

end