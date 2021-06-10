function display_sysInfo(sysInfo)
%DISPLAY_SYSINFO 此处显示有关此函数的摘要
fprintf('Settings of space mesh:      L= %2.1f, M=%i, dx = %2.3f\n', sysInfo.L, sysInfo.M, sysInfo.dx);
fprintf('Settings of time mesh:       T= %2.2f, dt = %2.3f\n', sysInfo.T, sysInfo.dt);
fprintf('Nonlinear function:          nlfn = %s \n', sysInfo.nlfn);
end

