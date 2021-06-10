function sysInfo = settings_DEs(display)
%%  Settings of the data

v = 1;                      % viscosity constant 
dt  = 0.001;                     % step size
T   = 1;                     % total time
L  = 10;                        % range of the basis
M = 300;                        % number of space grids

Initial = 'DN_1_0.5';
% Available Initial values:
% Initial = 'DN_2_0.5';
% Initial = 'N_0_1';
% Initial = 'DiffN_3_0.5_-1_1';
% Initial = 'DiffN_3_1_-1_1_0_1';
% Initial = 'unif_0.3';


nlfn = 'power_3';
% Available nonlinear functions:

% nlfn = 'Opinion_P1_3_-1_P2_4_2';
% nlfn = 'Opinion_P1_1_2_P2_3_-1_P3_4_2';
% nlfn = 'Opinion_P1_1_2_P2_3_-0.3_P3_4_0.1'
% nlfn = 'Opinion_P1_1_2_P2_2_0.5'
% nlfn = 'Opinion_P1_0.5_3_P2_2_-2_P3_2.5_0_P4_4_-5_P5_6_6';
% nlfn = 'Opinion_P1_1_3_P2_2_0.5_P3_3_-1_P4_4_-3_P5_6_6';

% nlfn = 'power_3';
% nlfn = 'LJ_3_0.2';
% nlfn = 'cos_5_2';
% nlfn = 'AB_2_0.5_threshold_0.05';
%% Save to structure 
sysInfo.v = v;
sysInfo.T = T;
sysInfo.dt = dt;
sysInfo.L = L;
sysInfo.M = M;
sysInfo.Initial = Initial;
sysInfo.nlfn = nlfn;

%% Get more information from the settings
sysInfo = get_sysInfo_details(sysInfo);
% 
% if nargin < 1
%     fprintf('Settings of space mesh:      L= %2.1f, M=%i, dx = %2.3f\n', L, M, sysInfo.dx);
%     fprintf('Settings of time mesh:       T= %2.2f, dt = %2.3f\n', T, dt);
%     fprintf('Nonlinear function:          nlfn = %s \n', nlfn);
% end

end
