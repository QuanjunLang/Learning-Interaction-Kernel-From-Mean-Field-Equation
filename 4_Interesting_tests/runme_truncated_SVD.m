% Givern the density u(x,t)
% Infer the interaction kernel phi

%% 
clear all;
close all;
clc
addtopath;

%% specify the settings for this script
sysInfo = settings_DEs();

sysInfo.dx = 0.025;
sysInfo.L = 10;
sysInfo.M = 400;
sysInfo.T = 0.5;
sysInfo.dt = 0.01;
sysInfo.Initial = 'DiffN_3_1_-1_1_0_1';
sysInfo.nlfn = '5cos(2x)';
sysInfo.v = 0.5;

sysInfo = get_sysInfo_details(sysInfo);

%% Generate/Load Data U
progressON = 0;
U = SPCC_Imp_solver(sysInfo, progressON);


%% Inference Computational kernel
% Computation kernel. Keep that, and the rest runs faster

G = inference_get_integration_kernel(sysInfo, U);

%% Inference (spline)
infer = settings_inference();
infer.type = 'spline_rho_support'; 
infer.b_method = 'ptgrad';
infer.c_method = 'Tik_L2_X_optimal_pinv';
infer.conv_method = 'fourier';
infer.rho_bdry = 5e-4;
infer.spline.deg = 0;
infer.spline.knot_num =  40;
infer.n = infer.spline.deg + infer.spline.knot_num;

infer_output = inference_kernel_all(infer, sysInfo, U, G);


Error = evaluation_error(sysInfo, infer_output, G);









%% Get some information for later use
A = infer_output.A;
b = infer_output.b;

bas = infer_output.basis.grad_K;

D_rho = zeros(infer.n, 1); % Discrete Rho
for i = 1:infer.n
    phi = bas(:, i);
    
    even_phi = abs(phi);
    D_rho(i) = trapz(trapz(U .* conv2(U, even_phi, 'same'))) * sysInfo.dt / sysInfo.T * sysInfo.dx;
end

mid = (sysInfo.M)/2+2;
rgrid = sysInfo.dx:sysInfo.dx:sysInfo.L;

D_rho_under_fine_grid = bas * D_rho;
D_rho_under_fine_grid = D_rho_under_fine_grid(mid:end);


temp = find(D_rho_under_fine_grid < (min(D_rho_under_fine_grid(D_rho_under_fine_grid>0))*2));
r_ind = temp(1);


ymin = min(infer_output.phi_kernel.true(mid:r_ind+mid))*1.1;
ymax = max(infer_output.phi_kernel.true(mid:r_ind+mid))*1.1;
yrange = ymax - ymin;
ymin = ymin - 0.3 * yrange;
ymax = ymax + 0.3 * yrange;


all_est = infer_output.phi_kernel.est;
% all_est = bas * (A\b.est);

FU = abs(fft(U)); FU = [FU(mid:end, :);FU(1:mid-1, :)];
% figure; surf(FU); shading flat;
%% RKHS with kernel G
[V, D] = eig(A);
[G_eigval, ind] = sort(diag(D),'descend');
G_eigvect = V(:, ind);





G_eigval_cumsum = cumsum(G_eigval);
G_eigval_cumsum = G_eigval_cumsum / G_eigval_cumsum(end);
all_cutoff = [0.9, 0.99, 0.999, 0.9999];
N = length(all_cutoff);

figure;
for i = 1:N
    cutoff = all_cutoff(i);
    
    ind = sum(G_eigval_cumsum < cutoff);
    
    subplot(2,N,i);
    plot(G_eigvect(:, 1:ind));
    title(['G cutoff = ',num2str(cutoff), ', ind = ', num2str(ind)],'FontSize',16)
    xlabel('radius r')
    set(gca,'FontSize',15);
    
    K = ind;
    AA = diag(G_eigval(1:K));
    bb = G_eigvect(:, 1:K)' * b.est;
    cc = G_eigvect(:, 1:K) * (AA\bb);
    
    G_est = bas * cc;
    
    
    subplot(2,N,i+N);
    
    yyaxis left 
    
    plot(rgrid, G_est(mid:end),'linewidth', 2);hold on;
    plot(rgrid, infer_output.phi_kernel.true(mid:end), 'r','linewidth',2)
    plot(rgrid, all_est(mid:end),'linewidth',2)
    xlim([sysInfo.dx, rgrid(r_ind)])
    ylim([ymin, ymax])
    ylabel('\phi')
    xlabel('radius r')
    
    yyaxis right
    fill(rgrid([1 1:end end]),[0 D_rho_under_fine_grid' 0],'k','EdgeColor',...
    'k','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2);
    ylabel('\rho')
    
    set(gca,'FontSize',15);
    
    legend('RKHS Estimation','True','all spline estimate', '\rho')
end
sgtitle('RKHS H_G Inference result','FontSize',20)
set(gcf,'unit','centimeters','position',[10 5 50 30]);

saveas(gcf,'G_infer.png')
%% RKHS with kernel R
P = diag(D_rho);

% [V, D] = eig(P\A);
[V, D] = eig(A, P);
[R_eigval, ind] = sort(diag(D),'descend');
R_eigvect = V(:, ind);



coef = diag(R_eigvect' * P * R_eigvect);
R_eigvect = R_eigvect ./ sqrt(coef');

R_eigval_cumsum = cumsum(R_eigval);
R_eigval_cumsum = R_eigval_cumsum / R_eigval_cumsum(end);
all_cutoff = [0.9, 0.99, 0.999, 0.9999];
N = length(all_cutoff);

figure;
for i = 1:N
    cutoff = all_cutoff(i);
    
    ind = sum(R_eigval_cumsum < cutoff);
    
    subplot(2,N,i);
    plot(R_eigvect(:, 1:ind));
    title(['R cutoff = ',num2str(cutoff), ', ind = ', num2str(ind)],'FontSize',16)
    xlabel('radius r')
    set(gca,'FontSize',15);
    
    K = ind;
    AA = R_eigvect(:, 1:K)' * A * R_eigvect(:, 1:K);
    bb = R_eigvect(:, 1:K)' * (b.est);
    cc = R_eigvect(:, 1:K) * (AA\bb);
    
    R_est = bas * cc;

    subplot(2,N,i+N);
    
    yyaxis left 
    
    plot(rgrid, R_est(mid:end),'linewidth', 2);hold on;
    plot(rgrid, infer_output.phi_kernel.true(mid:end), 'r','linewidth',2)
    plot(rgrid, all_est(mid:end),'linewidth',2)
    xlim([sysInfo.dx, rgrid(r_ind)])
    ylim([ymin, ymax])
    ylabel('\phi')
    xlabel('radius r')
    
    yyaxis right
    fill(rgrid([1 1:end end]),[0 D_rho_under_fine_grid' 0],'k','EdgeColor',...
    'k','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2);
    ylabel('\rho')
    
    set(gca,'FontSize',15);
    
    legend('Estimation','True','all spline', '\rho')
end
sgtitle('RKHS H_R Inference result','FontSize',20)
set(gcf,'unit','centimeters','position',[10 5 50 30]);
saveas(gcf,'R_infer.png')




%% b decay analysis
figure;
subplot(221)
plot(log10([abs(G_eigval), abs(G_eigvect' * b.est)]));legend('\lambda','\beta^G');
ylim([-30, 10])
title('log10 plot of the decay speed of \lambda and \beta^G')
subplot(222)
plot(log10([abs(R_eigval), abs(R_eigvect' *  b.est)]));legend('\gamma','\beta^R')
ylim([-30, 10])
title('log10 plot of the decay speed of \gamma and \beta^R')

subplot(223)
plot(log10(abs(G_eigvect' * b.est))./log10(abs(G_eigval)));title('relative decaying speed of G (\beta^G_i = \lambda_i^p)')
ylabel('p')


subplot(224)
plot(log10(abs(R_eigvect' *  b.est)) ./ log10(abs(R_eigval)));title('relative decaying speed of R (\beta^R_i = \gamma_i^p)')
ylabel('p')

sgtitle('Relative decay compare of \beta and eigenvalues')
saveas(gcf,'b_decay.png')


figure; 
semilogy(abs(G_eigval),'k-.x','linewidth', 1); hold on; 
semilogy(abs(R_eigval),'b-x','linewidth', 1); hold on; 
semilogy(abs(G_eigvect' * b.est),'k-.o','linewidth', 2); hold on; 
semilogy(abs(R_eigvect' * b.est),'b-.o','linewidth', 2); 
legend('G eigenvalues', 'R eigenvalues','G: b-projection', 'R: b-projection')

h= figure;  % compare spectrum vs b-prorjection decay
subplot(211);
semilogy(abs(G_eigval),'k-.x','linewidth', 1); hold on; 
semilogy(abs(R_eigval),'b-x','linewidth', 1); hold on; 
semilogy(abs(G_eigvect' * b.est),'k-.','linewidth', 1.5); hold on; 
semilogy(abs(R_eigvect' * b.est),'b-','linewidth', 1.5); 
legend('G eigenvalues', 'R eigenvalues','G: b-projection', 'R: b-projection')

subplot(212)
% figure %  decay of ratio  b-proj / spectrum
semilogy(abs(G_eigvect' * b.est)./abs(G_eigval),'k-.','linewidth',1); hold on;  
semilogy(abs(R_eigvect' * b.est)./abs(R_eigval),'b-','linewidth',1);
legend('Ratio: G','Ratio: R','location','best'); title('Ratio  b-projection/spectrum');
xlabel('Eigenvalue number'); 
set_positionFontsAll; 
saveas(h,'Decay_of_ratio.pdf'); 

%% b Error analysis
e = b.est-b.best;
figure;
yyaxis left; plot(abs(e));hold on;yyaxis right;plot(D_rho)
saveas(gcf,'b_error.png')

figure;
subplot(221)
plot(log10([abs(G_eigval), abs(G_eigvect' * e)]));legend('eigval','b');
ylim([-30, 10])
title('G')
subplot(222)
plot(log10([abs(R_eigval), abs(R_eigvect' *  e )]));legend('eigval','b')
ylim([-30, 10])
title('R')
sgtitle('error in b decay speed compare')

subplot(223)
plot(log10(abs(G_eigvect' * e))./log10(abs(G_eigval)));title('relative decaying speed of G')

subplot(224)
plot(log10(abs(R_eigvect' *  e)) ./ log10(abs(R_eigval)));title('relative decaying speed of G')

saveas(gcf,'b_error_decay.png')