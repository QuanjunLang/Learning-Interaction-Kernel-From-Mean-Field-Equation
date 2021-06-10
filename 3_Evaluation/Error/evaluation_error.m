function Error = evaluation_error(sysInfo, infer_output, G)
% get all kinds of error 
%%
dx = sysInfo.dx;
xgrid = sysInfo.xgrid;

rho = infer_output.rho;
phi_kernel = infer_output.phi_kernel;
A = infer_output.A;
b = infer_output.b;
c = infer_output.c;
B = infer_output.B;
lambda = infer_output.lambda;

%% L2_rho
Error.L2rho.best = sqrt(trapz(((phi_kernel.best-phi_kernel.true).^2.*rho.val_fine))*dx);
Error.L2rho.est = sqrt(trapz(((phi_kernel.est-phi_kernel.true).^2.*rho.val_fine))*dx);

%% L2_X
X = (abs(xgrid) < rho.support)';

Error.L2X.best = sqrt(trapz(((phi_kernel.best-phi_kernel.true).^2.*X))*dx);
Error.L2X.est = sqrt(trapz(((phi_kernel.est-phi_kernel.true).^2.*X))*dx);

%% E: use the true (not possible in practice)
Error.E.orical.best = (phi_kernel.best-phi_kernel.true)' * G * (phi_kernel.best-phi_kernel.true) * dx * dx;
Error.E.orical.est = (phi_kernel.est-phi_kernel.true)' * G * (phi_kernel.est-phi_kernel.true) * dx * dx;

%% E: in practice
Error.E.practical.best =  c.best' * A * c.best - 2* b.best' * c.best;
Error.E.practical.est =  c.est' * A * c.est - 2* b.est' * c.est;

%% E regularized
Error.E.regularized.best = c.best' * (A + lambda * B) * c.best - 2* b.best' * c.best;
Error.E.regularized.est = c.est' * (A + lambda * B) * c.est - 2* b.est' * c.est;
end