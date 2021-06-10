function b = inference_get_b(basis, sysInfo, infer, U, G)
% Generate the Nonhomogeneous term b

%% Extract parameters
basis_grad_K = basis.grad_K;
basis_K = basis.K;
basis_Hess_K = basis.Hess_K;

n = infer.n;
b_method = infer.b_method;
conv_method = infer.conv_method;

grad_K_true = sysInfo.phi_kernel;
L  = sysInfo.L;
dx = sysInfo.dx;
dt = sysInfo.dt;
v  = sysInfo.v;

xgrid = -L:dx:L;
%% 
b.best = zeros(n,1);
b.est  = zeros(n,1);




R_grad = function_convolution_with_U(basis_grad_K, U, dx, conv_method);
R_potential = function_convolution_with_U(basis_K, U, dx, conv_method);

%   b.best
for i = 1:n
    b.best(i) = grad_K_true(xgrid)*G*basis_grad_K(:,i)*dx*dx;
    
end

%   different methods for b
switch b_method
    case'ptgrad'
        for i = 1:n
            helper = trapz(gradient(U,dt).*R_potential(:,:,i))*dx;
            temp_b = -helper - v * trapz(gradient(U',dx)'.*R_grad(:,:,i))*dx;   % viscos
            b.est(i) = trapz(temp_b)*dt;
        end
    case'ptlap'
        R_laplace = function_convolution_with_U(basis_Hess_K, U, dx, conv_method);
        for i = 1:n
            helper = trapz(gradient(U,dt).*R_potential(:,:,i))*dx;
            temp_b = -helper + v * trapz(U.*R_laplace(:,:,i))*dx;     % viscos
            b.est(i) = trapz(temp_b)*dt;
        end
end

end