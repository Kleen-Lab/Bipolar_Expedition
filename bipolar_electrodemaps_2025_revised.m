
% Script to plot stacked histograms for number of bipolar pairs
% stratified by electrode type


data_root = getenv("BIPOLAR_DATA");
data_results = getenv("RESULTS");
%Grids/strips in blue, depths in red


% Create stacked histograms of elec distances by type

figure; set(gcf,'color','w','position',[100 100 1200 1000]);

pts = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', ...
    'EC219', 'EC220', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC162', 'EC168'}; % 16 patients used in analysis


skip_grids = {'EC181', 'EC220', 'EC162'}; % no subdural components
skip_depths = {};
skip_strips = {'EC181', 'EC220', 'EC162'}; % no subdural components


colors = [0.12 0.47 0.71;
    0.60 0.20 0.35;
    0.99 0.75 0.18;
    0.58 0.40 0.74;
    0.55 0.76 0.29;
    0.30 0.69 0.83;
    0.80 0.36 0.36;
    0.96 0.51 0.75;
    0.50 0.50 0.50;
    0.18 0.55 0.55;
    0.42 0.68 0.48;
    0.84 0.71 0.42;
    0.69 0.45 0.62;
    0.47 0.75 0.68;
    0.94 0.60 0.32;
    0.62 0.49 0.78];

electrode_types = {'grids', 'depths', 'strips'};
titles = {'Grid Electrodes', 'Depth Electrodes', 'Strip Electrodes'};
skip_lists = {skip_grids, skip_depths, skip_strips};

bin_edges = 0:1:80;
bin_centers = 0.5:1:79.5;

for e = 1:3

    subplot('Position', [0.1, 0.72 - (e-1)*0.3, 0.65, 0.25]);
    all_counts = zeros(length(bin_centers), length(pts));

    for p = 1:length(pts)

        if ismember(pts{p}, skip_lists{e})
            continue;
        end

        filename = [fullfile(data_results,['fix_' electrode_types{e} '_dist_tent_' pts{p} '.mat'])];

        if exist(filename, 'file')
            load(filename, 'distance');
            dists = distance(~isnan(distance) & distance > 0);
            counts = histcounts(dists, bin_edges);
            all_counts(:, p) = counts';
        end
    end

    h = bar(bin_centers, all_counts, 'stacked', 'EdgeColor', 'none', 'BarWidth', 1);

    for p = 1:length(pts)
        h(p).FaceColor = colors(p, :);
    end

    xlabel('Distance (mm)');
    ylabel('Counts');
    title(titles{e});
    xlim([0 80]);
    grid on;

end

legend_entries = {};

for p = 1:length(pts)
    legend_entries{p} = ['Pt. ' num2str(p)];
end

lgd = legend(h, legend_entries, 'Position', [0.78, 0.3, 0.15, 0.4]);

