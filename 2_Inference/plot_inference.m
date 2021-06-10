function plot_inference(sysInfo, infer, infer_output, G, all_dir, saveON)
% generate the figures of spline and rkhs results
%% parameters
rho = infer_output.rho.val_fine;
dx = sysInfo.dx;

%% relative error, spline basiss
Error_spline = evaluation_error(sysInfo, infer_output, G);

phi_kernel = infer_output.phi_kernel;

L2_norm = sqrt(trapz(phi_kernel.true.^2 .* rho)* dx);
L2_err  = Error_spline.L2rho.est;


RKHS_norm = sqrt(phi_kernel.true' * G * phi_kernel.true * dx * dx);
RKHS_err  = sqrt(Error_spline.E.orical.est);

fprintf('RKHS norm = %2.2f \n', RKHS_norm);
fprintf('L2 rho norm = %2.2f \n', L2_norm);

fprintf('RKHS relative error = %2.2f%% \n', RKHS_err/RKHS_norm*100)
fprintf('L2 rho relative error = %2.2f%% \n', L2_err/L2_norm*100)


%% -------- plotting -------

xgrid = sysInfo.xgrid;
mid = sysInfo.mid;
LW = 4;
L = sysInfo.L;

threshold = 0.01;
ind_temp = find(rho > threshold);
ind_temp = ind_temp(end);
rho_supp = xgrid(ind_temp);

fig = figure;
left_color = [0 0.4470 0.7410];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);


yyaxis right
fill([dx,dx:dx:L,L],[0 rho(mid+1:end)' 0],'k','EdgeColor','k','FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2);ylabel('\rho(r)')

% phi
yyaxis left
plot(dx:dx:L, infer_output.phi_kernel.true(mid+1:end),'-', 'Linewidth',LW);hold on;
plot(dx:dx:L, infer_output.phi_kernel.est(mid+1:end),'--','color', [0.8500 0.3250 0.0980], 'Linewidth',LW);
ylabel('\phi(r)'); xlabel('space radius r'); legend('True','Spline','\rho','location','best');xlim([dx rho_supp])

set(gca,'FontSize',10 );
fig = gcf;
fig.Units = 'inches';
fig.Position = [2 2   12 10];
% 
figname = [all_dir.fig_dir, sysInfo.nlfn, '_phi'];
setfigureFonts_saveFig(change_dot_to_p(figname), saveON)
end