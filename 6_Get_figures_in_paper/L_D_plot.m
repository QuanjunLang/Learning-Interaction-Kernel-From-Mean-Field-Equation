function L_D_plot(all_dir)
%%
load get_lambda_plot.mat
load get_dim_plot.mat
% 
%% Lambda L curve
X = [zeta, ita];
[~,R2,K2] = curvature(X);
sgn = -1 + 2 * (K2(:,1)>0).*(K2(:,2)>0);


[val, ind] = max(1./R2.*sgn);
lambda = all_lambda(ind);

figure;h_fig = subplot(121);
h = plot(zeta,ita, '-', 'LineWidth',4); grid on; axis equal;  set(h,'marker','.');
ylabel log_{10}(c'Bc);  xlabel log_{10}(c'Ac - 2b'c);    title('2D curve with curvature vectors');   hold on
quiver(zeta,ita,K2(:,1),K2(:,2),'LineWidth',1.5);   hold off
title('L-curve and normal vector');
legend('L-curve','normal vector')
sgn = -1 + 2 * (K2(:,1)>0).*(K2(:,2)>0);
set(h_fig, 'position', [0.13 0.1 0.3 0.8] )
set(gca,'FontSize',25 );
ylim([1,4.5])
xlim([-3.5,0])

h_fig = subplot(122);
plot(temp(ind), val, '*','Color', '#D95319','MarkerSize',15, 'LineWidth',5 );hold on;
plot(temp, 1./R2 .*sgn,'Color','#0072BD', 'LineWidth',4);
xlabel('log_{10}(\lambda)'); ylabel('curvature');legend(['\lambda_0 = 10^{', sprintf('%1.2f', temp(ind)), '}'],'Location','northeast')
title('signed curvature')
xlim([min(temp), max(temp)])


set(h_fig, 'position', [0.56 0.23 0.32 0.535] )


fig = gcf;
fig.Units = 'inches';
fig.Position = [2 2   14 8];

set(gca,'FontSize',25 );
figname = [all_dir.fig_dir,'L'];
print([figname,'.eps'],'-depsc');
%% Dimension

figure(11); % plot costFn, and optimal dimension
plot(all_basis_num,Err_real_est, ':o', 'LineWidth',4);hold on;
plot(all_basis_num,Err.real_reg.est, '--x', 'LineWidth',4);
title('The cost function E and E_\lambda');hold on;
xlabel('Dimension n');ylabel('Cost function');
plot(optimal_basis_num, val_EL, '*', 'MarkerSize',20, 'LineWidth',5);legend('E','E_\lambda',['Optimal dimension = ' ,num2str(optimal_basis_num)] )
set(gca,'FontSize', 15);
fig = gcf; fig.Units = 'inches';     fig.Position = [2 2   10 6];
set(gca,'FontSize',25 );
figname = [all_dir.fig_dir,'D'];
print([figname,'.eps'],'-depsc');
end