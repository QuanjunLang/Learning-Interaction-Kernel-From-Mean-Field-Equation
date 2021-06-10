function free_Engy = evaluation_Free_energy(sysInfo, infer_output, U)
% Compute the free energy along the solution

%% Load parameters
dx = sysInfo.dx;
L = sysInfo.L;

mid = sysInfo.mid;
xgrid = -L:dx:L;

phi_kernel = infer_output.phi_kernel;
%% interaction potential: est and true
Phi_potential.true = sysInfo.Phi_potential(xgrid');
Phi_potential.true = Phi_potential.true - Phi_potential.true(mid);
temp = cumsum(phi_kernel.est)*dx;
Phi_potential.est = temp - temp(mid);

%% Free energy along the data u
free_Engy.est  = free_energy(U, Phi_potential.est, dx, sysInfo.v);
free_Engy.true = free_energy(U, Phi_potential.true, dx, sysInfo.v);


end