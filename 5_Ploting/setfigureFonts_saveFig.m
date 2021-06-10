%% set figure fonts,
function setfigureFonts_saveFig(figname,saveON, WLratio,yrange,lgnd,location)
% setfigureFonts_saveFig(lgnd{1:nLagRun},'best',K,[2,1,1],figpath,figname)
%{
INPUT:
    figname   - figname without extension;
    WLratio   - width/Length ratio; [4,3,1] or [2,1,1];
    yrange    - ylim range
    lgnd,location - legend and its location.
%}
if exist('lgnd','var');  legend(lgnd,'Location',location); end
% xlim([0.8,K+0.5]); xticks(1:K);
if exist('yrange','var'); ylim(yrange); end

ftsz = 40; mksz = 30;
set(findall(gcf,'-property','FontSize'),'FontSize',ftsz);
% set(findall(gcf,'-property','MarkerSize'),'MarkerSize',mksz);
if exist('WLratio','var'); pbaspect(WLratio); end
%    yticklabels({'-1.0','0.0','1.0','2.0','3.0'}); yticklabels({'-2.0','-1.0','0.0','1.0','2.0'});
% savefig([figpath,figname,'.fig']);
if saveON
    print([figname,'.eps'],'-depsc');
end
end