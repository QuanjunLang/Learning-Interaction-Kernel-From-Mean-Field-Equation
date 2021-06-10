function infer  = settings_inference()
% settings of inference
%% basis
basis_type = 'spline';
basis_num = 'deg_2_knotnum_26';
% spline:
%   'spline_rho_support'
%   'spline'

% rkhs: 
%   rkhs:             Reproducing Kernel Hilbert Space, using only odd basis.
%   rkhs_odd_cutoff:  RKHS, odd basis and curoff for the kernel
%   rkhs_rho:         Fei's new basis
%   rkhs_denormal:    RKHS without normalization using eigenvalues.
%   rkhs_rad:         Use Radical kernel to generate eigenvectors

%% methods
b_method = 'ptgrad'; 
% ptgrad:          use one step of integration by parts
% ptlap:           use two steps of integration by parts

c_method = 'Tik_H1_unif_optimal_pinv';
%Adb:               c = A\b
%AtA                c = (A'*A)\(A'*b)
%pinv:              c = pinv(A)*b

conv_method = 'fourier';
% conv_same:         conv2(f,U,'same')
% fourier:           fourier transform fft

rho_bdry = 5e-4; % the threhold for the support set of rho

%% Assemble
infer.type = basis_type;
infer.basis_num = basis_num;

infer.b_method = b_method;
infer.c_method = c_method;
infer.conv_method = conv_method;
infer.rho_bdry = rho_bdry;

infer = get_infer_details(infer);
end
  