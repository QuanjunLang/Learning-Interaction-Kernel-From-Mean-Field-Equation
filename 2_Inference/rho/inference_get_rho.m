function rho = inference_get_rho(sysInfo, U, threshold)
% Generate the rho corresponding to the data U
%% Extract parameters
M = sysInfo.M;
TN = sysInfo.TN;
dx = sysInfo.dx;
L = sysInfo.L;

threshold = threshold * dx;
%%
rho_t = zeros(2*M+1,TN+1);
for k = 1:M+1
    rho_t(k,:)       = sum(U(1:k,:).*U(M+2-k:M+1,:))*dx;
    rho_t(2*M+2-k,:) = rho_t(k,:);
end
rho_t = rho_t(M/2+1:3*M/2+1,:);
val   = mean(rho_t,2);


% %% 
% rho_t = zeros(M+1, TN+1);
% U_flip = fliplr(U')';
% for t = 1:TN+1
%     rho_t(:, t) = conv(U_flip(:, t), U(:, t), 'same')*dx;
% end

%% Support width for rho
ind     = find(val>threshold);
assert(length(ind)>1);
ind    = ind(end);
rho_support = (-L + ind * dx);

%% rho with respect to radius
% half = val(M/2+1:end);
% half = half*2;
% half(1) = half(1)/2;
%% export
rho.val_fine = val;
rho.support = rho_support;
rho.support_index = ind;
% rho.val_fine_half = half;
end