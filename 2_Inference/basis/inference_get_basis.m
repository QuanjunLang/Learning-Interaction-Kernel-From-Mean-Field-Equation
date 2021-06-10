function basis = inference_get_basis(sysInfo, infer, rho, U, G)
% generate basis for inference

%% Extract parameters
% dx = sysInfo.dx;
type = infer.type;
% M = sysInfo.M;
%% switch case

if contains(type, 'spline')
    basis = get_spline_basis(sysInfo, infer, rho);
elseif contains(type, 'rkhs')
    basis = get_rkhs_basis(sysInfo, infer, G, rho);
end




end















