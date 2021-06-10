function plot_U(U, sysInfo, all_dir, saveON)
% plot the data U

% %%
xgrid = sysInfo.xgrid;
tgrid = sysInfo.tgrid;
% 
% [XX,YY] = meshgrid(tgrid, xgrid);
% 
% figure;
% surf(XX,YY,U); shading flat;colorbar;  view(0,90);
% xlabel('time t');ylabel('space x'); 
% set_positionFontsAll; 
% 
% title('Data U')

%%
threshold = 1e-30;
ind_U = find(min(U')' > threshold);
ind_U = ind_U(end);
U_supp = xgrid(ind_U);

figure;[XX, YY] = meshgrid(tgrid, xgrid);
surf(XX, YY, U); shading flat;view(0, 90);colorbar;
xlabel('time t');ylabel('space x');ylim([-U_supp, U_supp])
% title('U');

fig = gcf;fig.Units = 'inches';
fig.Position = [2 2   12 10];

figname = [all_dir.fig_dir, sysInfo.nlfn, '_U'];

setfigureFonts_saveFig(change_dot_to_p(figname), saveON)
end

