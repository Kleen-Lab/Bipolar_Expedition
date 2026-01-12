
% Code to plot figure 4: relative changes from referential 
% for each vs. all analysis, using data saved (.mat files) and generated from
% EachVsAll_cleaned2025.m

binz_re = 0:2:84;
frx_re = 2:2:200;

pts_grids = {'EC133', 'EC175', 'EC183', 'EC186', 'EC187', 'EC196', ...
'EC219', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC168'};
pts_names = {'Pt. 1', 'Pt. 2', 'Pt. 4', 'Pt. 5', 'Pt. 6', ...
'Pt. 7', 'Pt. 8', 'Pt. 10', 'Pt. 11', 'Pt. 12', 'Pt. 13', ...
'Pt. 14', 'Pt. 16'};
subplot_pos = [2 3 4 6 7 8 10 11 12 13 14 15 16];

figure;

all_mDiff = nan(100, 43, length(pts_grids));

data_root = getenv("BIPOLAR_DATA");
data_results = getenv("RESULTS");


% Load all mDiff matrices
for i = 1:length(pts_grids)

    subplot(4,4,subplot_pos(i))

    load(fullfile(data_results,['fix_grids_tent_' pts_grids{i} '.mat']));
    
    % Store this patient's mDiff in the 3D array
    all_mDiff(:,:,i) = mDiff;
    
    pcolorjk(binz_re,frx_re,mDiff)
    colormap(cmocean('balance',64))
    caxis([-50 50]);
    ylim([2 80.001]);
    xlim([4 60]);
    set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
    set(gca,'xtick',[10 20 30 40 50 60],'xticklabel',{'10','20','30','40','50','60'})
    title(pts_names{i});
end

% Compute the nanmean across the 3rd dimension (across patients)
mDiff_mean = nanmean(all_mDiff, 3);

save(fullfile(data_results,'m_eva_grids.mat'), 'mDiff_mean');


figure;

pcolorjk(binz_re,frx_re,mDiff_mean)
colormap(cmocean('balance',128))
caxis([-50 50]);
ylim([2 80.001]);
xlim([4 20.001]);
set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
set(gca,'xtick',[4:2:20],'xticklabel',{'4','6','8','10','12','14','16','18','20'})
xlabel('Distance (mm)');
ylabel('Frequency (Hz)');
grid on;
set(gca, 'Layer', 'top');
set(gca, 'GridColor', 'k');        % Black grid lines
set(gca, 'GridAlpha', 0.5);        % 50% opacity (darker)
set(gca, 'MinorGridLineStyle', 'none'); 
title('MEAN GRIDS')
colorbar;

%% STRIPS

% Initialize a 3D array to store all mDiff matrices
%pts_strips = {'EC133', 'EC175', 'EC183', 'EC186', 'EC187', 'EC196', ...
%'EC219', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC168'};
%pts_names = {'Pt. 1', 'Pt. 2', 'Pt. 4', 'Pt. 5', 'Pt. 6', ...
%'Pt. 7', 'Pt. 8', 'Pt. 10', 'Pt. 11', 'Pt. 12', 'Pt. 13', ...
%'Pt. 14', 'Pt. 16'};
pts_strips = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', ...
'EC219', 'EC220', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC162', 'EC168'};
pts_names = {'Pt. 1', 'Pt. 2', 'Pt. 3', 'Pt. 4', 'Pt. 5', 'Pt. 6', ...
'Pt. 7', 'Pt. 8', 'Pt. 9', 'Pt. 10', 'Pt. 11', 'Pt. 12', 'Pt. 13', ...
'Pt. 14', 'Pt. 15', 'Pt. 16'};

all_mDiff_strips = nan(100, 22, length(pts_strips));

figure;
% Load all mDiff matrices from strips
for i = 1:length(pts_strips)
    
    subplot(4,4,i)
    title(pts_names{i});

    if ~(strcmp(pts_strips{i}, 'EC181') || strcmp(pts_strips{i}, 'EC220') || strcmp(pts_strips{i}, 'EC162')) 
        load(fullfile(data_results,['fix_strips_tent_' pts_strips{i} '.mat']));
        
        % Store this patient's mDiff in the 3D array
        all_mDiff_strips(:,:,i) = mDiff;
        
        pcolorjk(binz_re,frx_re,mDiff)
        colormap(cmocean('balance',64))
        caxis([-50 50]);
        ylim([2 80.001]);
        xlim([4 60]);
        set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
        set(gca,'xtick',[10 20 30 40 50 60],'xticklabel',{'10','20','30','40','50','60'})
        title(pts_names{i});
    end
end

% Compute the nanmean across the 3rd dimension (across patients)
mDiff_mean = nanmean(all_mDiff_strips, 3);

% Plot the mean
figure;
pcolorjk(binz_re,frx_re,mDiff_mean)
colormap(cmocean('balance',128))
caxis([-50 50]);
ylim([2 80.001]);
xlim([4 20.001]);
set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
set(gca,'xtick',[4:2:20],'xticklabel',{'4','6','8','10','12','14','16','18','20'})
xlabel('Distance (mm)');
ylabel('Frequency (Hz)');

% Add darker grid on top - major grid only
grid on;
set(gca, 'Layer', 'top');
set(gca, 'GridColor', 'k');
set(gca, 'GridAlpha', 0.5);
set(gca, 'MinorGridLineStyle', 'none');
title('MEAN STRIPS')

colorbar;

%% DEPTHS

% Initialize a 3D array to store all mDiff matrices
pts_depths = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', ...
'EC219', 'EC220', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC162', 'EC168'};
pts_names = {'Pt. 1', 'Pt. 2', 'Pt. 3', 'Pt. 4', 'Pt. 5', 'Pt. 6', ...
'Pt. 7', 'Pt. 8', 'Pt. 9', 'Pt. 10', 'Pt. 11', 'Pt. 12', 'Pt. 13', ...
'Pt. 14', 'Pt. 15', 'Pt. 16'};

all_mDiff_depths = nan(100, 43, length(pts_depths));
figure;
% Load all mDiff matrices from depths
for i = 1:length(pts_depths)
    subplot(4,4,i)
    load(fullfile(data_results,['fix_depths_tent_' pts_depths{i} '.mat']));
    
    % Store this patient's mDiff in the 3D array
    all_mDiff_depths(:,:,i) = mDiff;
    
    pcolorjk(binz_re,frx_re,mDiff)
    colormap(cmocean('balance',64))
    caxis([-50 50]);
    ylim([2 80.001]);
    xlim([4 60]);
    set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
    set(gca,'xtick',[10 20 30 40 50 60],'xticklabel',{'10','20','30','40','50','60'})
    title(pts_names{i});
end

% Compute the nanmean across the 3rd dimension (across patients)
mDiff_mean = nanmean(all_mDiff_depths, 3);

% Plot the mean
figure;
pcolorjk(binz_re,frx_re,mDiff_mean)
colormap(cmocean('balance',128))
caxis([-50 50]);
ylim([2 80.001]);
xlim([4 20.001]);
set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
set(gca,'xtick',[4:2:20],'xticklabel',{'4','6','8','10','12','14','16','18','20'})
xlabel('Distance (mm)');
ylabel('Frequency (Hz)');

% Add darker grid on top - major grid only
grid on;
set(gca, 'Layer', 'top');
set(gca, 'GridColor', 'k');
set(gca, 'GridAlpha', 0.5);
set(gca, 'MinorGridLineStyle', 'none');

title('MEAN DEPTHS')
colorbar;

%% BIPOLAR SQRT

figure;
pcolorjk(binz_re,frx_re,mARb_m);
colormap('jet');
ylim([2 80.001]);
xlim([4 60]);
set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
set(gca,'xtick',[10 20 30 40 50 60],'xticklabel',{'10','20','30','40','50','60'})
title('Referential');

figure;
pcolorjk(binz_re,frx_re,mb_m);
colormap('jet');
ylim([2 80.001]);
xlim([4 60]);
set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
set(gca,'xtick',[10 20 30 40 50 60],'xticklabel',{'10','20','30','40','50','60'});
title('Bipolar');




%%

binz_re = 0:2:84;
frx_re = 2:2:200;

% GRIDS
pts_grids = {'EC133', 'EC175', 'EC183', 'EC186', 'EC187', 'EC196', ...
'EC219', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC168'};

all_mDiff_grids = nan(100, 43, length(pts_grids));

for i = 1:length(pts_grids)
    load(fullfile(data_results,['fix_grids_tent_' pts_grids{i} '.mat']));
    all_mDiff_grids(:,:,i) = mDiff;
end

mDiff_mean_grids = nanmean(all_mDiff_grids, 3);

% STRIPS
pts_strips = {'EC133', 'EC175', 'EC183', 'EC186', 'EC187', 'EC196', ...
'EC219', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC168'};

all_mDiff_strips = nan(100, 43, length(pts_strips));

for i = 1:length(pts_strips)
    load(fullfile(data_results,['fix_depths_tent_' pts_strips{i} '.mat']));
    all_mDiff_strips(:,:,i) = mDiff;
end

mDiff_mean_strips = nanmean(all_mDiff_strips, 3);

% DEPTHS
pts_depths = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', ...
'EC219', 'EC220', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC162', 'EC168'};

all_mDiff_depths = nan(100, 43, length(pts_depths));

for i = 1:length(pts_depths)
    load(fullfile(data_results,['fix_depths_tent_' pts_depths{i} '.mat']));
    all_mDiff_depths(:,:,i) = mDiff;
end

mDiff_mean_depths = nanmean(all_mDiff_depths, 3);

% Create 1x3 subplot
figure('Position', [100 100 1400 400]);

% Plot 1: MEAN GRIDS
subplot(1,3,1)
pcolorjk(binz_re, frx_re, mDiff_mean_grids)
colormap(cmocean('balance',128))
caxis([-50 50]);
ylim([2 80.001]);
xlim([4 20.001]);
set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
set(gca,'xtick',[4:2:20],'xticklabel',{'4','6','8','10','12','14','16','18','20'})
xlabel('Distance (mm)');
ylabel('Frequency (Hz)');
grid on;
set(gca, 'Layer', 'top');
set(gca, 'GridColor', 'k');
set(gca, 'GridAlpha', 0.5);
set(gca, 'MinorGridLineStyle', 'none'); 
title('Grids (mean)', 'FontSize',16)
colorbar;

% Plot 2: MEAN STRIPS
subplot(1,3,3)
pcolorjk(binz_re, frx_re, mDiff_mean_strips)
colormap(cmocean('balance',128))
caxis([-50 50]);
ylim([2 80.001]);
xlim([4 20.001]);
set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
set(gca,'xtick',[4:2:20],'xticklabel',{'4','6','8','10','12','14','16','18','20'})
xlabel('Distance (mm)');
ylabel('Frequency (Hz)');
grid on;
set(gca, 'Layer', 'top');
set(gca, 'GridColor', 'k');
set(gca, 'GridAlpha', 0.5);
set(gca, 'MinorGridLineStyle', 'none');
title('Strips (mean)', 'FontSize',16)

% Plot 3: MEAN DEPTHS
subplot(1,3,2)
pcolorjk(binz_re, frx_re, mDiff_mean_depths)
colormap(cmocean('balance',128))
caxis([-50 50]);
ylim([2 80.001]);
xlim([4 20.001]);
set(gca,'yscale','log','ytick',[2 5 10 20 40 80],'yticklabel',{'2','5','10','20','40','80'})
set(gca,'xtick',[4:2:20],'xticklabel',{'4','6','8','10','12','14','16','18','20'})
xlabel('Distance (mm)');
ylabel('Frequency (Hz)');
grid on;
set(gca, 'Layer', 'top');
set(gca, 'GridColor', 'k');
set(gca, 'GridAlpha', 0.5);
set(gca, 'MinorGridLineStyle', 'none');
title('Depths (mean)', 'FontSize',16)



