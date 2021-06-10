function sysInfo = get_sysInfo_details(sysInfo)
% Get more information from simple sysInfo settings
%%
dt = sysInfo.dt;
T = sysInfo.T;
L = sysInfo.L;
M = sysInfo.M;

%%
tgrid = 0:dt:T;             % time grid

TN       = length(tgrid)-1;    % time grid number
dx = L * 2 /M;
xgrid = -L:dx:L;
rgrid = 0:dx:L;



U0 = get_Initial(L, dx, sysInfo.Initial);
[phi_kernel, Phi_potential] = get_nlfn(sysInfo.nlfn, L, dx);

%%
sysInfo.xgrid = xgrid;
sysInfo.tgrid = tgrid;
sysInfo.rgrid = rgrid;
sysInfo.mid = M/2+1;
sysInfo.TN = TN;
sysInfo.dx = dx;
sysInfo.U0 = U0;
sysInfo.Phi_potential = Phi_potential;
sysInfo.phi_kernel = phi_kernel;

end