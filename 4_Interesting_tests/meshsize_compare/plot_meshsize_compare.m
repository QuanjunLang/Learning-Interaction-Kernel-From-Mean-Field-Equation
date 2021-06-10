function plot_meshsize_compare(sysInfo, all_infer_output, all_M, all_dir, saveON)
%PLOT_MESHSIZE_COMPARE 此处显示有关此函数的摘要
%% Compute error
% use the finest mesh to compute error
N = length(all_M);
L = sysInfo.L;
finest_A = all_infer_output{end}.A;
finest_b = all_infer_output{end}.b.est;
finest_xgrid = -L:(2*L/all_M(end)):L;
finest_dx = 2*L/all_M(end);
finest_phi_true = sysInfo.phi_kernel(finest_xgrid);
finest_rho = all_infer_output{end}.rho.val_fine';

Error.E = zeros(N, 1);
Error.L2rho = zeros(N, 1);

for kM = 1:N
    small_xgrid = -L:(2*L/all_M(kM)):L;
    est = all_infer_output{kM}.phi_kernel.est;
    temp = fit(small_xgrid', est ,'linearinterp');
    est_func_handle = @(x) reshape(temp(x), size(x));
    c = all_infer_output{kM}.c.est;
    
    
    Error.L2rho(kM) = sqrt(trapz((est_func_handle(finest_xgrid) - finest_phi_true).^2.*finest_rho)*finest_dx);
    
    
    Error.E(kM) = c'*finest_A*c -2*c'*finest_b;
    
end


a = 1;







%% plot
% --------------- L2 rho error ------------
LW = 4;
all_dx = reshape(2*L./all_M, 1, []);
all_L2rho = reshape(Error.L2rho, 1, []);
switch sysInfo.nlfn
    case 'power_3'
        ind = 1:10;
        opt_L2rho_slope = 2;
    case 'AB_2_-0.5_cut_threshold_0.05'
        ind = [1,2,4,5,6,7,8,9,10];
        opt_L2rho_slope = 4/3;
    case 'Opinion_P1_3_-1_P2_4_2'
        ind = [1,2,4,5,6,7,8,9,10];
        opt_L2rho_slope = 2;
end

all_dx = all_dx(ind);
all_L2rho = all_L2rho(ind);


p = polyfit(log10(all_dx),log10(all_L2rho),1);
yfit = polyval(p, log10(all_dx));
fun2 = @(dx) 10^p(2) * dx.^p(1);

figure;
loglog(all_dx, all_L2rho, 'o','Linewidth',LW, 'Markersize',10);hold on;
loglog(all_dx, fun2(all_dx),'Linewidth',LW);
loglog(all_dx(1), fun2(all_dx(1)), '.','Linewidth',0.001, 'Markersize', 0.001)

box on
grid on
legend('Test point',sprintf('Slope = %1.2f',p(1)),sprintf('Optimal = %1.2f',opt_L2rho_slope),'Location','best')
xlabel('\Delta x');  ylabel('$L^2(\bar\rho_T)$ error','Interpreter','Latex')


fig = gcf; fig.Units = 'inches'; fig.Position = [2 2   12 10];
figname = [all_dir.fig_dir, sysInfo.nlfn, '_L2rE'];
setfigureFonts_saveFig(change_dot_to_p(figname), saveON)

%% -------------- Error functional E using optimization ------
all_dx = reshape(2*L./all_M, 1, []);
all_E = reshape(Error.E, 1, []);
switch sysInfo.nlfn
    case 'power_3'
        ind = 1:10;
        opt_E_slope = 4;
    case 'AB_2_-0.5_cut_threshold_0.05'
        ind = [1,2,4,5,6,7,8,9,10];
        opt_E_slope = 8/3;
    case 'Opinion_P1_3_-1_P2_4_2'
%         ind = [1,3,4,5,6,7,8,9,10];
        ind = [1,2,3,4,5,6,7,8,9,10];
        opt_E_slope = 4;
end

all_dx = all_dx(ind);
all_E = all_E(ind);

% x0 = [4.9,6.5,5];
% cost2 = @(x) sum((all_E - (x(3)*all_dx.^x(1) - x(2))).^2);
% options = optimset('Display','iter','PlotFcns',@optimplotfval);
% x = fminsearch(cost2,x0);
% fun1 = @(dx) x(3)*dx.^x(1) - x(2);


switch sysInfo.nlfn
    case 'power_3'
        x0 = [5,6.9,5];
        lb = [0, -min(all_E)/2, 0];
        ub = [5,7, 10];
    case 'AB_2_-0.5_cut_threshold_0.05'
        x0 = [5,2.74,5];
        lb = [0, -min(all_E)/2, 0];
        ub = [5, 3, 10];
    case 'Opinion_P1_3_-1_P2_4_2'
        x0 = [1,0.435,5];
        lb = [0, -min(all_E)/2, 0];
        ub = [3,0.46, 10];
end

A = [];b = [];Aeq = [];beq = [];
cost1 = @(x) sum((log10(all_E + x(2)) - x(1)*log10(all_dx) - log10(x(3))).^2);

nonlcon = [];
options = optimoptions('fmincon','Display','off');

x = fmincon(cost1,x0,A,b,Aeq,beq,lb,ub, nonlcon, options);
fun1 = @(dx) x(3)*dx.^x(1) - x(2);
% figure;plot(all_dx, fun1(all_dx));hold on;scatter(all_dx, all_E)


figure;

loglog(all_dx, all_E + x(2), 'o','Linewidth',LW, 'Markersize',10);hold on;
loglog(all_dx, fun1(all_dx) + x(2),'Linewidth',LW);
loglog(all_dx(1), fun1(all_dx(1)) + x(2),'.','Markersize',0.001,'Linewidth',0.001)

temp = gca;old_Ytick = temp.YTick;
% ------- using the original number, with certain digits of numbers ---
% 
% new_Ytick = split(sprintf('%1.4f ',old_Ytick - x(2)));
% yticklabels(new_Ytick);

% ------- using c + 10^e ---

new_Ytick = cell(length(old_Ytick));
for i = 1:length(old_Ytick)
    new_Ytick{i} = sprintf('%1.2f + 10^{%d}',-x(2), log10(old_Ytick(i)));
end
yticklabels(new_Ytick);

% -------------------------

legend('Test point',sprintf('Slope = %1.2f',x(1)),sprintf('Optimal = %1.2f',opt_E_slope),'Location','best')
xlabel('\Delta x');ylabel('Error functional E')
grid on


fig = gcf; fig.Units = 'inches'; fig.Position = [2 2   12 10];
figname = [all_dir.fig_dir, sysInfo.nlfn, '_EE_opt'];
setfigureFonts_saveFig(change_dot_to_p(figname), saveON);
% 
% figure
% scatter(all_dx, all_E, 90,'Linewidth',LW);hold on;
% plot(all_dx, fun1(all_dx),'Linewidth',LW);





end

