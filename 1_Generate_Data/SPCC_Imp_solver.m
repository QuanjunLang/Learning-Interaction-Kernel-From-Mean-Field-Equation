function U = SPCC_Imp_solver(sysInfo, progressON)
%Use explicit SPCC to solve Nonlinear Fokker-Planck equation
%% Extract the parameters
dx = sysInfo.dx;
dt = sysInfo.dt;
M = sysInfo.M;
TN = sysInfo.TN;
P = sysInfo.phi_kernel;
U0 = sysInfo.U0;
v = sysInfo.v;


if progressON
    fs={'\b\b%d%%','\b\b\b%d%%'};%waitbar
    fprintf('Generate Data U, SPCC Implicit Solver, progress:0%%');
end
%% Main
%Convolution helper function
C_mat = Convolution_Matrix(M, P, dx);

%Initialize U
U = zeros(M+1, TN+1);
f0 = U0;
U(:,1) = f0;
for i = 1: TN
    
    lambda = trapz(C_mat.*f0)'*dx/v;                      %(M+2)*1, [-1,0,1,2,....M]
    delta = 1./lambda + 1./(1-exp(lambda));             %(M+2)*1, [-1,0,1,2,....M]
    delta(delta==inf)=1/2;
    
    C = lambda/dx*v;                                      %(M+2)*1, [-1,0,1,2,....M]
    A = dt/dx*(C.*delta-v/dx);                          %(M+2)*1, [-1,0,1,2,....M]
    B = dt/dx*(C.*(1-delta)+v/dx);                      %(M+2)*1, [-1,0,1,2,....M]
    
    Mat = [zeros(1,M+1);diag(A(2:end-1)),zeros(M,1)];
    Mat = Mat + eye(M+1);
    Mat = Mat - diag([A(2:end-1);0]);
    Mat = Mat + diag([0;B(2:end-1)]);
    Mat = Mat + [zeros(M,1), diag(-B(2:end-1));zeros(1,M+1)];
    
    f0 = Mat\f0;
    U(:,i+1) = f0;
    if isnan(sum(f0))
        fprintf("SPCC solution blows up at time %f\n", i*dt)
        break
    end
    if progressON
        j = floor(i/(TN/100));
        fprintf(fs{1+(j>=10)},j);
    end
end
if progressON; fprintf('\n'); end

end

function C_mat = Convolution_Matrix(M,P,dx)
    function i = ind(k)
        i = k + M + 2;
    end
C = zeros(2*M+2, 1);


for i = -M-1:1:M
    C(ind(i)) = integral(P,i*dx, i*dx+dx);
end


C_mat = zeros(M+1,M+2);
for i = 1:M+1
    C_mat(i,:) = C(ind(-i:M+1-i));
end

C_mat(isinf(C_mat)) = P(dx)*dx;
end

