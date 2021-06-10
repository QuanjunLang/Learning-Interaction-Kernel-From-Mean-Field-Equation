function G = inference_get_integration_kernel(sysInfo, U)
% get basis for RKHS
dx = sysInfo.dx;
dt = sysInfo.dt;

[xsize, tsize] = size(U);
M = xsize - 1;
TN = tsize - 1;
L = M*dx/2;
T = TN*dt;

function temp = U_ext(ind)
    assert((ind >=-M/2)&&(ind <= M/2));
    if ind == 0
        temp = U;
    elseif ind > 0
        temp = [U(ind+1:end,:);zeros(ind, TN+1)];
    else
        temp = [zeros(-ind, TN+1);U(1:end+ind,:)];
    end
end

G = zeros(M+1,M+1);
fprintf('Progress: '); reverseStr = []; 
for i = 1:M+1
    reverseStr = displayprogress(100*i/(M+1), reverseStr);
    for j = 1:i
        G(i,j) = trapz(trapz(U.*U_ext(i-1-M/2).*U_ext(j-1-M/2)))*dx*dt;
        G(j,i) = G(i,j);
    end
end

end








