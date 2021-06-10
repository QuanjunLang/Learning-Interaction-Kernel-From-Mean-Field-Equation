function A = inference_get_A(basis, sysInfo, infer, G)
% Generate the normal matrix A
%%
n = infer.n;
dx = sysInfo.dx;
basis_grad_K = basis.grad_K;

A = zeros(n,n);
for i = 1:n
    for j = 1:i
        A(i,j) = basis_grad_K(:, i)' * G * basis_grad_K(:,j) *dx*dx;
        A(j,i) = A(i,j);
    end
end
end