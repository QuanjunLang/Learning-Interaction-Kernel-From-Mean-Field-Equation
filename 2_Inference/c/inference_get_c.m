function [phi_kernel, c] = inference_get_c(basis, sysInfo, infer, A, B, b, lambda)
% From A, B, lambda, b get c
% c = (A + lambda B)\b
% and return the estimated grad K

%% Extract parameters
basis_val_grad_K = basis.grad_K;
c_method = infer.c_method;
grad_K_true = sysInfo.phi_kernel;
L  = sysInfo.L;
dx = sysInfo.dx;
xgrid = -L:dx:L;
%%
if contains(c_method, 'Tik')
    reg_string = strsplit(c_method,'_');
    devide_method = reg_string{5};
else
    devide_method = c_method;
end


AB = A + lambda * B;


switch devide_method
    case 'Adb'
        c.est = AB\b.est;
        c.best = A\b.best;
    case'AtA'
        c.est = (AB'*AB)\(AB'*b.est);
        c.best = (A'*A)\(A'*b.best);
    case 'pinv'
        c.est = pinv(AB)*b.est;
        c.best = pinv(A)*b.best;
end

phi_kernel.true = grad_K_true(xgrid)';
phi_kernel.est = basis_val_grad_K*c.est;
phi_kernel.best = basis_val_grad_K*c.best;

phi_kernel.est = (phi_kernel.est - fliplr(phi_kernel.est')')/2;
phi_kernel.best = (phi_kernel.best - fliplr(phi_kernel.best')')/2;
end