
function fig3_plot

%SCRIPT TO PLOT FIG 3
% this script plots two patients using ECXXX_fig3_data generated from the
% fig3_EachVsAll.m plot

data_dir = getenv("BIPOLAR_DATA");
data_results = getenv("RESULTS");

sizeoffont = 12;
xldist = [0 60];
yldist = [2 80];
frx = 2:2:200;
xt = [10 20 30 40 50 60]; xtl = cellstr(num2str(xt'));

ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft'));

% make sure to add helper function paths


fig = figure(5);
set(fig, 'Position', [100, 100, 1400, 900]);
set(gcf, 'Color', 'white');

load(fullfile(data_results,'EC175_fig3_data.mat'));

subplot1 = subplot('Position', [0.375, 0.68, 0.24, 0.24]);
%elecsbrain(pt,0,[1:256],[0 0 0],'l',0,3.3,2); alpha 1; % requires montages

subplot2 = subplot('Position', [0.37, 0.25, 0.25, 0.4]);
pcolorjk(binz(1:size(toplot,2)),frx,toplot); shading flat; set(gca,'ydir','normal'); 
ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',sizeoffont); %colorbar;
%title({pt ' ln(power), z-scored by frequency',''},'fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl, 'xtick', xt, 'xticklabel', xtl);
xlim(xldist); ylim(yldist);
clim([-1 1]*(max(abs(clim)))); colormap(gca,cmocean('curl')); %cm=cmocean('balance',100); colormap(gca,cm(:,[2 1 3]))
clim([-4.1 4.1]);

subplot3 = subplot('Position',[0.37, 0.13, 0.25, 0.06]);
histogram(make1d(Mbp_distance),[0.001 binsz:binsz:85],'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); 
%xlabel('Binned bipolar distance (mm)','fontsize',sizeoffont);
ylabel('# pairs','fontweight','normal'); 
axis tight; grid on; %cb=colorbar; set(cb,'visible','off'); 
xlim(xldist); set(gca,'fontsize',sizeoffont);

clear pt binz toplot frx binsz Mbp_distance cm_distance

% EC183 ADDING ON

load(fullfile(data_results,'EC183_fig3_data.mat'));

subplot4 = subplot('Position', [0.67, 0.68, 0.25, 0.25]);
%elecsbrain(pt,0,[1:256],[0 0 0],'l',0,3.3,2); alpha 1; % requires montages

subplot5 = subplot('Position', [0.67, 0.25, 0.284, 0.4]);
pcolorjk(binz(1:size(toplot,2)),frx,toplot); shading flat; set(gca,'ydir','normal'); 
ylabel('Frequency (Hz)'); xlabel('Distance (mm)'); set(gca,'fontsize',sizeoffont); %colorbar;
%title({pt ' ln(power), z-scored by frequency',''},'fontweight','normal')
set(gca,'yscale','log','ytick',ft,'yticklabel',ftl, 'xtick', xt, 'xticklabel', xtl);
xlim(xldist); ylim(yldist)
clim([-1 1]*(max(abs(clim)))); colormap(gca,cmocean('curl')); %cm=cmocean('balance',100); colormap(gca,cm(:,[2 1 3]))
clim([-4.1 4.1]); cb = colorbar; set(cb, 'visible', 'on');

subplot6 = subplot('Position',[0.67, 0.13, 0.25, 0.06]);
histogram(make1d(Mbp_distance),[0.001 binsz:binsz:85],'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); 
%xlabel('Binned bipolar distance (mm)','fontsize',sizeoffont);
ylabel('# pairs','fontweight','normal'); 
axis tight; grid on; %cb=colorbar; set(cb, ); 
xlim(xldist); set(gca,'fontsize',sizeoffont);

end



