function plot_Free_energy(free_Engy, sysInfo, all_dir, saveON)

dt = sysInfo.dt;
T = sysInfo.T;
LW = 4;
%% --------- Free energy ---------
figure;
free_Engy.est_const = free_Engy.est - (free_Engy.est(1) - free_Engy.true(1));
plot(0:dt:T, free_Engy.true,':','Linewidth',LW*2); hold on;
plot(0:dt:T, free_Engy.est_const,'--','Linewidth',LW*2);
xlabel('Time t'); ylabel('Free energy')
legend('True', 'Estimated')
box on
% title('Free Energy of U(\cdot, t)')
fig = gcf; fig.Units = 'inches'; fig.Position = [2 2   12 10];
figname = [all_dir.fig_dir, sysInfo.nlfn, '_FE'];
setfigureFonts_saveFig(change_dot_to_p(figname), saveON)
end