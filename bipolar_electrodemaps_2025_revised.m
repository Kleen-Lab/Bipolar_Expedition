
% Script to plot patient brains and electrodes (Figure 1) for 
% the 16 patients in the analysis.
% Additionally plot stacked histograms for number of bipolar pairs
% stratified by electrode type

function bipolar_electrodemaps_2025_revised

    fig = figure(1);
    set(fig, 'Position', [100, 100, 1400, 900]);
    set(gcf, 'Color', 'white');
    pts = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', ...
    'EC219', 'EC220', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC162', 'EC168'};
    pts_names = {'Pt. 1', 'Pt. 2', 'Pt. 3', 'Pt. 4', 'Pt. 5', 'Pt. 6', ...
    'Pt. 7', 'Pt. 8', 'Pt. 9', 'Pt. 10', 'Pt. 11', 'Pt. 12', 'Pt. 13', ...
    'Pt. 14', 'Pt. 15', 'Pt. 16'};
    
    % Depth elec numbers pulled directly from TDT info sheet.
    
    depths = {[281:320],
     [289:340],
     [],
     [257:286 289:298],
     [257:266 289:318],
     [299:318],
     [273:282 289:308],
     [399:408 429:468],
     [],
     [353:382 385:404],
     [257:298],
     [271:300 309:328], 
     [321:344],
     [321:340], 
     [1:30 33:62 65:94 97:116], 
     [327:346 353:372]}; 
    
    %Manually specifying to create  tighter subplots
    
    h_spacing = 0.02; 
    v_spacing = 0.03; 
    h_margin = 0.03;  
    v_margin = 0.04;  
    
    n_rows = 4;
    n_cols = 4;
    subplot_width = (1 - 2*h_margin - (n_cols-1)*h_spacing) / n_cols;
    subplot_height = (1 - 2*v_margin - (n_rows-1)*v_spacing) / n_rows;
    
    %Grids/strips in blue, depths in red
    
    for i = 1:length(pts)
        disp(pts{i})
        
        row = ceil(i/n_cols);
        col = mod(i-1, n_cols) + 1;
        
        left = h_margin + (col-1)*(subplot_width + h_spacing);
        bottom = 1 - v_margin - row*subplot_height - (row-1)*v_spacing;
        
        ax = subplot('Position', [left, bottom, subplot_width, subplot_height]);
        
    
        if strcmp(pts{i}, 'EC220')
            elecsbrain(pts{i}, 0, [], [0 0 1], 'b', 0, 2.75, 2); lightsout; litebrain('i', 1); alpha 0.3;
            elecsbrain(pts{i}, 0, depths{i}, [1 0 0], 'b', 0, 3.8, 2); lightsout; litebrain('i', 1); alpha 0.3;
            view(-180,-90);
        else
            elecsbrain(pts{i}, 0, [], [0 0 1], 'b', 0, 2.75, 2); lightsout; litebrain('l', 1); alpha 0.3;
            elecsbrain(pts{i}, 0, depths{i}, [1 0 0], 'b', 0, 3.8, 2); lightsout; litebrain('l', 1); alpha 0.3;
        end
        
        title(pts_names{i});
    end
    
    %
    
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
        
            filename = ['/data' ...
            'fix_' electrode_types{e} '_dist_tent_' pts{p} '.mat'];
            
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

end
