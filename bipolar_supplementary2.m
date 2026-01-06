
% PLOT SUPPLEMENTARY FIGURE 2
% Display changes in IED spatial and temporal extent for HD vs
% subsampled conditions using different linelength window sizes

% Script pulls directly from data saved from high_density_ecog.m script

output_chs = [];
output_widths = [];

conds = {'LL20', 'LL40', 'LL100', 'absDer'};

for c = 1:size(conds,2)

    % path to LL20.mat, LL40.mat, LL100.mat, absDer.mat 
    % NOTE: these files are created when running
    % high_density_ecog_script.m

    load(['/data/results/' conds{c} '.mat']);

    strx = permResultsCell;
    num_chs = zeros(1,size(strx,2));
    mean_width = zeros(1,size(strx,2));
    
    
    for i = 1:size(strx,2)
        hd_chs = strx{i}.numSigChannels;
        subsamp_chs = strx{i}.numSigChannelsSub;
        hd_width = strx{i}.meanWidth;
        subsamp_width = strx{i}.meanWidthSub;
        num_chs(i) = hd_chs - subsamp_chs;
        mean_width(i) = hd_width - subsamp_width;
    end

    output_chs = [output_chs; num_chs];
    output_widths = [output_widths; mean_width];

    disp(['DONE: ' conds(c)]);

end

mean_width_new = output_widths;
num_chs_new = output_chs;


figure;
condition_labels = {'LL20', 'LL40', 'LL100', 'Absolute Derivative'};

subplot(2,1,1);
violinplot_with_lines(mean_width',condition_labels)


title('Differences in IED Duration', 'FontSize', 14);
ylabel('mean Width (ms)', 'FontSize',13);
set(gcf, 'Color', 'w');

subplot(2,1,2);
violinplot_with_lines(num_chs',condition_labels)


title('Differences in IED Spatial Extent', 'FontSize', 14);
ylabel('# of Channels Involved', 'FontSize',13);
set(gcf, 'Color', 'w');

