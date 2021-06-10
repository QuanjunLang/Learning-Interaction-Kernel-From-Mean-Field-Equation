function plot_Wasserstein_distance(ws_dist, sysInfo, all_dir, saveON)



%%
dt = sysInfo.dt;
T = sysInfo.T;
%% plotting
figure;
LW = 4;
plot(0:dt:T, [ws_dist.old_initial,ws_dist.new_initial],'Linewidth',LW)
% title('Wasserstein d(U_{true}(\cdot,t), U_{est}(\cdot,t))')
xlabel('Time t')
ylabel('Wasserstein distance')
legend('Original initial','New initial')
fig = gcf;
fig.Units = 'inches';
fig.Position = [2 2   12 10];

figname = [all_dir.fig_dir, sysInfo.nlfn, '_Wd'];
setfigureFonts_saveFig(change_dot_to_p(figname), saveON)

end