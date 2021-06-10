function [B, lambda] = inference_get_regularization(basis, sysInfo, infer, U, rho, A, b, plot_L_curve)
% Use A and b to get c. Choose normalization methods.
basis_val_grad_K = basis.grad_K;
basis_val_Hess_K = basis.Hess_K;
n = infer.n;
c_method = infer.c_method;
L  = sysInfo.L;
dx = sysInfo.dx;
xgrid = -L:dx:L;
%%
if contains(c_method, 'Tik')
    reg_string = strsplit(c_method,'_');
    % Measure of the H norm
    switch reg_string{3}
        case 'unif'
            density = ones(size(xgrid));
        case 'rho'
            density = rho.val_fine;
        case 'X'
            density = xgrid < rho.support;
    end
    
    % H0 or H1 norm
    switch reg_string{2}
        case 'H0'
            B = basis_val_Hess_K' * diag(density) * basis_val_Hess_K *dx;
        case 'H1'
            B = basis_val_Hess_K' * diag(density) * basis_val_Hess_K *dx;
            B = B + basis_val_grad_K' * diag(density) * basis_val_grad_K;
        case 'I'
            B = eye(n);
        case 'L2'
            B = basis_val_grad_K' * diag(density) * basis_val_grad_K;
    end
    
    % find regularization parameter lambda
    if strcmp(reg_string{4}, 'optimal')
        lambda = get_optimal_lambda(A, B, b.est, plot_L_curve);
    else
        lambda = str2double(reg_string{4});
    end
    
else
    B = eye(n);
    lambda = 0;
end








end