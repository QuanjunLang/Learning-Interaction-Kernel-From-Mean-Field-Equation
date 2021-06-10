function basis = get_rkhs_basis(sysInfo, infer, G, rho)
% Get RKHS basis

%% Load parameters
n = infer.n;
temp = split(infer.type, '+');
type = temp{1};
dx = sysInfo.dx;
M = sysInfo.M;
mid = sysInfo.mid;
switch type
    case{'rkhs'}
        [V,~] = svd(G);
        basis.grad_K = V(:, 1:n);
    case{'rkhs_rho'}
        [V, D] = eig(G, diag(rho.val_fine));
        [~, ind] = sort(diag(D), 'descend');
        V = V(:, ind);
        basis.grad_K = V(:, 1:n);
    case{'rkhs_mesh'}
        temp_infer = infer;
        temp_infer.type = 'spline';
        temp_infer.basis_num = temp{2};
        temp_infer = get_infer_details(temp_infer);
        
        temp_basis = get_spline_basis(sysInfo, temp_infer, rho);
        
        A = inference_get_A(temp_basis, sysInfo, temp_infer, G);
        [V, ~] = svd(A);
        basis.grad_K = temp_basis.grad_K * V(:, 1:n);
    case{'rkhs_rho_mesh'}
        temp_infer = infer;
        temp_infer.type = 'spline';
        temp_infer.basis_num = temp{2};
        temp_infer = get_infer_details(temp_infer);
        
        temp_basis = get_spline_basis(sysInfo, temp_infer, rho);
        
        A = inference_get_A(temp_basis, sysInfo, temp_infer, G);
        P = diag(-rho.val_fine(1:mid)' * temp_basis.grad_K(1:mid, :));
        [V, D] = eig(A, P);
        [~, ind] = sort(diag(D), 'descend');
        V = V(:, ind);
        
        basis.grad_K = temp_basis.grad_K * V(:, 1:n);
end


basis.Hess_K = gradient(basis.grad_K',dx)';
basis.K = cumsum(basis.grad_K)*dx;
basis.K = basis.K - basis.K(M/2+1,:);

end