function infer_output = inference_kernel_all(infer, sysInfo, U, G, plot_L_curve)
% Inference the interaction kernel from data U

% density rho for the pairwise difference (regression measure)
rho_threshold = 0.001;
rho = inference_get_rho(sysInfo, U, rho_threshold);

% Get basis information
basis = inference_get_basis(sysInfo, infer, rho, U, G);

% Get b
b = inference_get_b(basis, sysInfo, infer, U, G);

% Get matrix A
A = inference_get_A(basis, sysInfo, infer, G);

% Get regularization matrix B and lambda
% plot_L_curve = 1;
[B, lambda] = inference_get_regularization(basis, sysInfo, infer, U, rho, A, b, plot_L_curve);

% Get c and grad_K
[phi_kernel, c] = inference_get_c(basis, sysInfo, infer, A, B, b, lambda);


%% return
infer_output.A = A;
infer_output.b = b;
infer_output.c = c;
infer_output.rho = rho;
infer_output.B = B;
infer_output.lambda = lambda;
infer_output.phi_kernel = phi_kernel;
infer_output.basis = basis;
end








%% Best possible result of the true grad_k, using given basis
% [grad_K_val.truth_approx, c.truth_approx]  = RKHS_get_best_func_by_density(grad_K_true, basis_val_grad_K, dx, L, rho); % best possible based on the given basis
% b.truth_approx = A*c.truth_approx;
% %% get grad_K
% if contains(c_method, 'Tik')
%     reg_string = strsplit(c_method,'_');
%     % Measure of the H norm
%     if strcmp(reg_string{3}, 'unif')
%         density = ones(size(xgrid));
%     elseif strcmp(reg_string{3}, 'rho')
%         density = rho;
%     elseif strcmp(reg_string{3}, 'inv')
%         density = U(:,end);
%     elseif strcmp(reg_string{3}, 'antirho')
%         density = 1 - rho;
%         density = density - min(density);
%     end
%     
%     % H0 or H1 norm
%     switch reg_string{2}
%         case 'H0'
%             B = basis_val_Hess_K' * diag(density) * basis_val_Hess_K *dx;
%         case 'H1'
%             B = basis_val_Hess_K' * diag(density) * basis_val_Hess_K *dx;
%             B = B + basis_val_grad_K' * diag(density) * basis_val_grad_K;
%         case 'I'
%             B = eye(n);
%     end
%     
%     % Load regularization parameter lambda
%     if strcmp(reg_string{4}, 'optimal')
%         %%%%%%%%% print? 1, 0
%         lambda = get_optimal_lambda(A, B, b.est, 0);
%     else
%         lambda = str2double(reg_string{4});
%     end
%     
%     
%     % Adb method
%     devide_method = reg_string{5};
%     
%     
%     [grad_K_val.best, c.best] = get_func(A,b.best,basis_val_grad_K, devide_method);                % using grad K true
%     [grad_K_val.est,c.est] = get_func(A,b.est,basis_val_grad_K, devide_method, B, lambda);
%     infer_output.B = B;
%     infer_output.reg = lambda;
%     
%     
%     % Make it odd
%     grad_K_val.best = (grad_K_val.best - fliplr(grad_K_val.best')')/2;
%     grad_K_val.est = (grad_K_val.est - fliplr(grad_K_val.est')')/2;
%     
%     
%     Error_E_real_reg.best = c.best' * (A + lambda * B) * c.best - 2* b.best' * c.best;
%     Error_E_real_reg.truth_approx = c.truth_approx' * (A + lambda * B) * c.truth_approx - 2* b.truth_approx' * c.truth_approx;
%     Error_E_real_reg.est = c.est' * (A + lambda * B) * c.est - 2* b.est' * c.est;
%     
%     
%     infer_output.Error_E_real_reg = Error_E_real_reg;
% else
%     [grad_K_val.est,c.est] = get_func(A,b.est,basis_val_grad_K, c_method);
%     [grad_K_val.best, c.best] = get_func(A,b.best,basis_val_grad_K, c_method);
% end
% 
% 
% grad_K_val.true = grad_K_true(xgrid)';
% 






% %% Error
% Error_L2rho.best = sqrt(sum(((grad_K_val.best'-grad_K_true(xgrid)).^2*rho))*dx);
% Error_L2rho.est = sqrt(sum(((grad_K_val.est'-grad_K_true(xgrid)).^2*rho))*dx);
% Error_L2rho.truth_approx = sqrt(sum(((grad_K_val.truth_approx'-grad_K_true(xgrid)).^2*rho))*dx);
% 
% Error_L2.best = sqrt(sum(((grad_K_val.best'-grad_K_true(xgrid)).^2))*dx);
% Error_L2.est = sqrt(sum(((grad_K_val.est'-grad_K_true(xgrid)).^2))*dx);
% Error_L2.truth_approx = sqrt(sum(((grad_K_val.truth_approx'-grad_K_true(xgrid)).^2))*dx);
% 
% 
% %fprintf('Error finished.\n')
% 
% 
% Conv.best = get_funcs_conv_U(grad_K_val.best, U, dx, R_method);
% Conv.est = get_funcs_conv_U(grad_K_val.est, U, dx, R_method);
% Conv.truth_approx = get_funcs_conv_U(grad_K_val.truth_approx, U, dx, R_method);
% Conv.True = get_funcs_conv_U(grad_K_true(xgrid'), U, dx, R_method);
% 
% 
% % Error_E.best            = trapz(trapz((Conv.best-Conv.True).^2.*U))*dx*dt;
% % Error_E.est             = trapz(trapz((Conv.est-Conv.True).^2.*U))*dx*dt;
% % Error_E.truth_approx    = trapz(trapz((Conv.truth_approx-Conv.True).^2.*U))*dx*dt;
% 
% Error_E.best            = (grad_K_val.best - grad_K_val.true)' * G * (grad_K_val.best - grad_K_val.true) * dx * dx;
% Error_E.est             = (grad_K_val.est - grad_K_val.true)' * G * (grad_K_val.est - grad_K_val.true) * dx * dx;
% Error_E.truth_approx    = (grad_K_val.truth_approx - grad_K_val.true)' * G * (grad_K_val.truth_approx - grad_K_val.true) * dx * dx;
% 


% 
% infer_output.Conv = Conv;
% infer_output.Error_E = Error_E;
% infer_output.err = Error_L2rho;
% infer_output.Error_L2 = Error_L2;
% %% Error for paper analysis
% inference_error = sqrt(sum(((grad_K_val.est'-grad_K_val.truth_approx').^2*rho))*dx);
% infer_output.inference_error = inference_error;
% %% Error_E_real
% Error_E_real.best = c.best' * A * c.best - 2* b.best' * c.best;
% Error_E_real.truth_approx = c.truth_approx' * A * c.truth_approx - 2* b.truth_approx' * c.truth_approx;
% Error_E_real.est = c.est' * A * c.est - 2* b.est' * c.est;
% 
% Error_E_real_plus_const.best = c.best' * A * c.best - 2* b.best' * c.best + grad_K_val.true' * G * grad_K_val.true *dx*dx;
% Error_E_real_plus_const.est = c.est' * A * c.est - 2* b.est' * c.est+ grad_K_val.true' * G * grad_K_val.true *dx*dx;
% Error_E_real_plus_const.truth_approx = c.truth_approx' * A * c.truth_approx - 2* b.truth_approx' * c.truth_approx+ grad_K_val.true' * G * grad_K_val.true *dx*dx;
% 
% 
% Reg_norm.best = c.best' * B * c.best *dx *dx;
% Reg_norm.est = c.est' * B * c.est *dx *dx;
% Reg_norm.truth_approx = c.truth_approx' * B * c.truth_approx *dx *dx;
% 
% 
% 
% infer_output.Reg_norm = Reg_norm;
% infer_output.Error_E_real = Error_E_real;
% infer_output.Error_E_real_plus_const = Error_E_real_plus_const;
% %% Save to output structure
% infer_output.L_adapt = L_adapt;
% 
% infer_output.rho = rho;
% 
% grad_K_val.basis = basis_val_grad_K;
% 
% infer_output.b = b;
% infer_output.A = A;
% infer_output.c = c;
% 
% basis.phi = basis_val_phi;
% basis.grad_K = basis_val_grad_K;
% basis.K = basis_val_K;
% basis.Hess_K = basis_val_Hess_K;
% 
% infer_output.basis = basis;
% infer_output.G = G;
% 
% infer_output.Conv = Conv;
% infer_output.Error_E = Error_E;
% 
% 
% infer_output.grad_K_val = grad_K_val;

% 
% %% function handel for possible cases
% if isa(basis_grad_K,'double')
%     infer_output.grad_K = grad_K_val;
% else
%     grad_K.truth_approx = c_times_basis(c.truth_approx, basis_grad_K);
%     grad_K.best = c_times_basis(c.best, basis_grad_K);
%     grad_K.est = c_times_basis(c.est, basis_grad_K);
%     grad_K.basis = basis_grad_K;
%     
%     infer_output.grad_K = grad_K;
% end



% end


% %% from A, b get the approximated grad K
% 
% function [k,c] = get_func(A, b, basis, c_method, B, lambda)
% 
% % regularzation term
% if nargin == 6
%     A = A + lambda * B;
% end
% 
% switch c_method
%     case 'Adb'
%         c = A\b;
%     case'AtA'
%         c = (A'*A)\(A'*b);
%     case 'pinv'
%         c = pinv(A)*b;
% end
% k = basis*c;
% end
% 
% 
% 
% 
% 
% 
% 
% 
% function temp = c_times_basis(c, basis)
% temp = @(x) 0*x;
% for i = 1:length(c)
%     temp = @(x) temp(x) + c(i) * basis{i}(x);
% end
% end




