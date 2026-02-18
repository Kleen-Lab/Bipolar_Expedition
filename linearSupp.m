
% Plot supplementary linear supplementary analysis for patients

data_dir = getenv("BIPOLAR_DATA");
result_dir = getenv("RESULTS");

% order of patients
pts = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', ...
       'EC219', 'EC220', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC162', 'EC168'};

start_pt = 1; % change to 1 for pts 1-4, 5 for pts 5-8, 9 for pts 9-12 etc.
end_pt = 4; %change to 4 for pts 1-4, 8 for pts 5-8, 12 for pts 9-12 etc, 


cm=cool(17);
cm=[0 0 0;1 1 1;1 1 1;cm];
frx = 2:2:200;
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft'));

% distances for each elec component type
depthsdist = [0, 5, 10, 15, 20];
gridsdist = [0, 4, 8, 12, 16, 20];
stripsdist = [0, 10, 20];
frxrange = [2 200];


% INSET plots
% ORDER: grid depth strip (columns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); set(gcf,'color','w','position',[372 1 1297 1337]);

comps = {'grid', 'depth', 'strip'};
trm_data_all = cell(1,3);
for c = 1:3
    load(fullfile(result_dir,[comps{c} '_trm_data.mat']));
    trm_data_all{c} = trm_data;
end

valid_pts = {};
for p = 1:length(pts)
    has_data = false;
    for c = 1:3
        pts_in_data = trm_data_all{c}.patients;
        pt_idx = find(strcmp(pts_in_data, pts{p}));
        if c == 1
            data = trm_data_all{1}.grids;
        elseif c == 2
            data = trm_data_all{2}.depths;
        else
            data = trm_data_all{3}.strips;
        end
        if ~isempty(pt_idx) && ~isempty(data{1,pt_idx})
            has_data = true;
            break;
        end
    end
    if has_data
        valid_pts{end+1} = pts{p};
    end
end

pts_to_plot = valid_pts(start_pt:min(end_pt, length(valid_pts)));

for i = 1:length(pts_to_plot)
    for c = 1:3 % 1=grids, 2=depths, 3=strips
        subplot_idx = (i-1)*3 + c;
        ax = subplot(length(pts_to_plot), 3, subplot_idx);
       
        if c == 1 % grids
            data = trm_data_all{1}.grids;
            dists = gridsdist;
            comp_name = 'Grid Electrodes';
        elseif c == 2 % depths
            data = trm_data_all{2}.depths;
            dists = depthsdist;
            comp_name = 'Depth Electrodes';
        else % strips
            data = trm_data_all{3}.strips;
            dists = stripsdist;
            comp_name = 'Strip Electrodes';
        end
        

        pts_in_data = trm_data_all{c}.patients;
        pt_idx = find(strcmp(pts_in_data, pts_to_plot{i}));
        
        if ~isempty(pt_idx) && ~isempty(data{1,pt_idx})
            % Plot ALL ribbons
            for j = 1:length(dists)

                if ~isempty(data{j,pt_idx}) && size(data{j,pt_idx}, 1) > 0 && size(data{j,pt_idx}, 2) > 0
                    ribbons(frx,data{j,pt_idx},cm(max([1 dists(j)]),:),.5,'sem',0,0);
                    set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl);
                    hold on;
                end
            end
            
            if c == 1 % grids
                legend({'referential','', [num2str(gridsdist(2)) ' mm'],'', [num2str(gridsdist(3)) ' mm'],'', [num2str(gridsdist(4)) ' mm'],'', [num2str(gridsdist(5)) ' mm'],'', [num2str(gridsdist(6)) ' mm']''}, 'location','sw');
            elseif c == 2 % depths
                legend({'referential','', [num2str(depthsdist(2)) ' mm'],'', [num2str(depthsdist(3)) ' mm'],'', [num2str(depthsdist(4)) ' mm'],'', [num2str(depthsdist(5)) ' mm']''}, 'location','sw');
            else % strips
                legend({'referential','', [num2str(stripsdist(2)) ' mm'],'', [num2str(stripsdist(3)) ' mm']''}, 'location','sw');
            end
            

            title(comp_name, 'FontWeight', 'normal', 'FontSize', 10);
            
            axis tight;
            set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl);
            grid on;
        else
      
            text(0.5,0.5,'No data','HorizontalAlignment','center');
            set(gca,'xlim',[0 1],'ylim',[0 1],'XTick',[],'YTick',[]);
            title(comp_name, 'FontWeight', 'normal', 'FontSize', 9);
        end
        

        if c == 2
            pos = get(ax, 'Position');

            annotation('textbox', [pos(1), pos(2)+pos(4)+0.01, pos(3), 0.03], ...
                'String', ['Pt. ' num2str(start_pt + i - 1)], ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 13, ...
                'FontWeight', 'bold', ...
                'EdgeColor', 'none', ...
                'FitBoxToText', 'off');
        end
    end
    
    % Display patient name in terminal after completing each row
    disp(['Pt. ' num2str(start_pt + i - 1) ' = ' pts_to_plot{i}]);
end

pts = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', ...
       'EC219', 'EC220', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC162', 'EC168'};

cm=cool(17);
cm=[0 0 0;1 1 1;1 1 1;cm];
frx = 2:2:200;
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft'));

%% Define distance arrays for each electrode type
depthsdist = [0, 5, 10, 15, 20];
gridsdist = [0, 4, 8, 12, 16, 20];
stripsdist = [0, 10, 20];
frxrange = [2 200];

%% ORDER: grid depth strip (columns)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); set(gcf,'color','w','position',[372 1 1297 1337]);
comps = {'grids', 'depths', 'strips'};

figure(2); set(gcf,'color','w','position',[372 1 1297 1337]);
comps = {'grid', 'depth', 'strip'};

trm_data_all = cell(1,3);
for c = 1:3
    load(fullfile(result_dir,[comps{c} '_trm_data.mat']));
    trm_data_all{c} = trm_data;
end

valid_pts = {};
for p = 1:length(pts)
    has_data = false;
    for c = 1:3
        pts_in_data = trm_data_all{c}.patients;
        pt_idx = find(strcmp(pts_in_data, pts{p}));
        if c == 1
            data = trm_data_all{1}.grids;
        elseif c == 2
            data = trm_data_all{2}.depths;
        else
            data = trm_data_all{3}.strips;
        end
        if ~isempty(pt_idx) && ~isempty(data{1,pt_idx})
            has_data = true;
            break;
        end
    end
    if has_data
        valid_pts{end+1} = pts{p};
    end
end


pts_to_plot = valid_pts(start_pt:min(end_pt, length(valid_pts)));


n_rows = length(pts_to_plot);
n_cols = 3;
h_margin = 0.05; 
v_margin = 0.06;  
h_spacing = 0.03; 
v_spacing = 0.06; 

subplot_width = (1 - 2*h_margin - (n_cols-1)*h_spacing) / n_cols;
subplot_height = (1 - 2*v_margin - (n_rows-1)*v_spacing) / n_rows;

% Plot each patient (rows) x each electrode type (columns)
for i = 1:length(pts_to_plot)
    for c = 1:3 % 1=grids, 2=depths, 3=strips

        row = i;
        col = c;
        left = h_margin + (col-1)*(subplot_width + h_spacing);
        bottom = 1 - v_margin - row*subplot_height - (row-1)*v_spacing;
        
        ax = axes('Position', [left, bottom, subplot_width, subplot_height]);
        
        if c == 1 % grids
            data = trm_data_all{1}.grids;
            dists = gridsdist;
            comp_name = 'Grid Electrodes';
        elseif c == 2 % depths
            data = trm_data_all{2}.depths;
            dists = depthsdist;
            comp_name = 'Depth Electrodes';
        else % strips
            data = trm_data_all{3}.strips;
            dists = stripsdist;
            comp_name = 'Strip Electrodes';
        end
        

        pts_in_data = trm_data_all{c}.patients;
        pt_idx = find(strcmp(pts_in_data, pts_to_plot{i}));
        
        if ~isempty(pt_idx) && ~isempty(data{1,pt_idx})
            % Plot ALL ribbons
            for j = 1:length(dists)
                if ~isempty(data{j,pt_idx}) && size(data{j,pt_idx}, 1) > 0 && size(data{j,pt_idx}, 2) > 0
                    ribbons(frx,data{j,pt_idx},cm(max([1 dists(j)]),:),.5,'sem',0,0);
                    set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl);
                    hold on;
                    xlabel('Frequency (Hz)', 'FontSize',7.5);
                    ylabel('sqrt(Power)', 'FontSize',7.5);
                end
            end

            if i == 1
                if c == 1 % grids
                    lgd = legend({'referential','', [num2str(gridsdist(2)) ' mm'],'', [num2str(gridsdist(3)) ' mm'],'', [num2str(gridsdist(4)) ' mm'],'', [num2str(gridsdist(5)) ' mm'],'', [num2str(gridsdist(6)) ' mm']''}, 'location','sw');
                elseif c == 2 % depths
                    lgd = legend({'referential','', [num2str(depthsdist(2)) ' mm'],'', [num2str(depthsdist(3)) ' mm'],'', [num2str(depthsdist(4)) ' mm'],'', [num2str(depthsdist(5)) ' mm']''}, 'location','sw');
                else % strips
                    lgd = legend({'referential','', [num2str(stripsdist(2)) ' mm'],'', [num2str(stripsdist(3)) ' mm']''}, 'location','sw');
                end
                lgd.FontSize = 7;
            end
            
            title(comp_name, 'FontWeight', 'normal', 'FontSize', 9);
            set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl);
            grid on;
            
            x_50_norm = (log10(50) - log10(2)) / (log10(200) - log10(2));
            x_200_norm = (log10(200) - log10(2)) / (log10(200) - log10(2));
            

            inset_left = left + subplot_width * x_50_norm;
            inset_width = subplot_width * (x_200_norm - x_50_norm);
            inset_height = subplot_height * 0.28 * 2.5;  
            inset_bottom = bottom + subplot_height * 0.68 - (inset_height - subplot_height * 0.28);  % Adjust position
            
            ax_inset = axes('Position', [inset_left, inset_bottom, inset_width, inset_height]);
            

            for j = 1:length(dists)
                if ~isempty(data{j,pt_idx}) && size(data{j,pt_idx}, 1) > 0 && size(data{j,pt_idx}, 2) > 0
                    ribbons(frx,data{j,pt_idx},cm(max([1 dists(j)]),:),.5,'sem',0,0);
                    hold on;
                end
            end
            
            set(gca,'xlim',[50 200],'xscale','log','xtick',[50 100 200],'FontSize',7);
            grid on;
            box on;
            
        else

            text(0.5,0.5,'No data','HorizontalAlignment','center');
            set(gca,'xlim',[0 1],'ylim',[0 1],'XTick',[],'YTick',[]);
            title(comp_name, 'FontWeight', 'normal', 'FontSize', 9);
        end
        

        if c == 2
            annotation('textbox', [left, bottom+subplot_height+0.005, subplot_width, 0.025], ...
                'String', ['Pt. ' num2str(start_pt + i - 1)], ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 11, ...
                'FontWeight', 'bold', ...
                'EdgeColor', 'none', ...
                'FitBoxToText', 'off');
        end
    end
    
    disp(['Pt. ' num2str(start_pt + i - 1) ' = ' pts_to_plot{i}]);
end


%% AGGREGATE PLOTS: Mean across all patients for each electrode type
% USING RIBBONS

figure(3); set(gcf,'color','w','position',[372 1 1297 400]);

comps = {'grid', 'depth', 'strip'};
trm_data_all = cell(1,3);
for c = 1:3
    load(fullfile(result_dir,[comps{c} '_trm_data.mat']));
    trm_data_all{c} = trm_data;
end

valid_pts = {};
for p = 1:length(pts)
    has_data = false;
    for c = 1:3
        pts_in_data = trm_data_all{c}.patients;
        pt_idx = find(strcmp(pts_in_data, pts{p}));
        if c == 1
            data = trm_data_all{1}.grids;
        elseif c == 2
            data = trm_data_all{2}.depths;
        else
            data = trm_data_all{3}.strips;
        end
        if ~isempty(pt_idx) && ~isempty(data{1,pt_idx})
            has_data = true;
            break;
        end
    end
    if has_data
        valid_pts{end+1} = pts{p};
    end
end

for c = 1:3
    ax = subplot(1, 3, c);
    pos = get(ax, 'Position');
    

    if c == 1 % grids
        data = trm_data_all{1}.grids;
        dists = gridsdist;
        comp_name = 'Grid Electrodes';
        pts_in_data = trm_data_all{1}.patients;
    elseif c == 2 % depths
        data = trm_data_all{2}.depths;
        dists = depthsdist;
        comp_name = 'Depth Electrodes';
        pts_in_data = trm_data_all{2}.patients;
    else % strips
        data = trm_data_all{3}.strips;
        dists = stripsdist;
        comp_name = 'Strip Electrodes';
        pts_in_data = trm_data_all{3}.patients;
    end
    

    all_patient_means = cell(length(dists), 1);
    
    for j = 1:length(dists)
        patient_means = [];
        for p = 1:length(valid_pts)
            pt_idx = find(strcmp(pts_in_data, valid_pts{p}));
            if ~isempty(pt_idx) && j <= size(data,1) && ~isempty(data{j,pt_idx}) && ...
                    size(data{j,pt_idx}, 1) > 0 && size(data{j,pt_idx}, 2) > 0
                
                %Taking mean across channels (rows) for this patient
                patient_mean = nanmean(data{j,pt_idx}, 1); 
                patient_means = [patient_means; patient_mean];
            end
        end
        all_patient_means{j} = patient_means;
        
        if ~isempty(patient_means)
            ribbons(frx, patient_means, cm(max([1 dists(j)]),:), .3, 'sem', 0, 0);
            set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl);
            hold on;
        end
    end
    

    if c == 1 % grids
        legend({'referential','', [num2str(gridsdist(2)) ' mm'],'', [num2str(gridsdist(3)) ' mm'],'', [num2str(gridsdist(4)) ' mm'],'', [num2str(gridsdist(5)) ' mm'],'', [num2str(gridsdist(6)) ' mm']''}, 'location','sw');
    elseif c == 2 % depths
        legend({'referential','', [num2str(depthsdist(2)) ' mm'],'', [num2str(depthsdist(3)) ' mm'],'', [num2str(depthsdist(4)) ' mm'],'', [num2str(depthsdist(5)) ' mm']''}, 'location','sw');
    else % strips
        legend({'referential','', [num2str(stripsdist(2)) ' mm'],'', [num2str(stripsdist(3)) ' mm']''}, 'location','sw');
    end
    
    title(comp_name, 'FontWeight', 'normal', 'FontSize', 10);
    xlabel('Frequency (Hz)');
    ylabel('sqrt(Power)');
    set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl);
    grid on;
    
    %inset
    x_50_norm = (log10(50) - log10(2)) / (log10(200) - log10(2));
    x_200_norm = (log10(200) - log10(2)) / (log10(200) - log10(2));
    
    inset_left = pos(1) + pos(3) * x_50_norm;
    inset_width = pos(3) * (x_200_norm - x_50_norm);
    inset_height = pos(4) * 0.70;  
    inset_bottom = pos(2) + pos(4) * 0.25;  
    
    ax_inset = axes('Position', [inset_left, inset_bottom, inset_width, inset_height]);
    

    for j = 1:length(dists)
        patient_means = all_patient_means{j};
        if ~isempty(patient_means)
            ribbons(frx, patient_means, cm(max([1 dists(j)]),:), .3, 'sem', 0, 0);
            hold on;
        end
    end
    
    set(gca,'xlim',[50 200],'xscale','log','xtick',[50 100 200],'FontSize',7);
    grid on;
    box on;
end

sgtitle('Aggregated: Mean Across All Patients', 'FontWeight', 'bold', 'FontSize', 14);

%% NOW MEANS USING MEANS AS LINES NOT RIBBONS

figure(4); set(gcf,'color','w','position',[372 1 1297 400]);

comps = {'grid', 'depth', 'strip'};
trm_data_all = cell(1,3);
for c = 1:3
    load(fullfile(result_dir,[comps{c} '_trm_data.mat']));
    trm_data_all{c} = trm_data;
end

valid_pts = {};
for p = 1:length(pts)
    has_data = false;
    for c = 1:3
        pts_in_data = trm_data_all{c}.patients;
        pt_idx = find(strcmp(pts_in_data, pts{p}));
        if c == 1
            data = trm_data_all{1}.grids;
        elseif c == 2
            data = trm_data_all{2}.depths;
        else
            data = trm_data_all{3}.strips;
        end
        if ~isempty(pt_idx) && ~isempty(data{1,pt_idx})
            has_data = true;
            break;
        end
    end
    if has_data
        valid_pts{end+1} = pts{p};
    end
end

for c = 1:3
    ax = subplot(1, 3, c);
    pos = get(ax, 'Position');
    
    if c == 1 % grids
        data = trm_data_all{1}.grids;
        dists = gridsdist;
        comp_name = 'Grid Electrodes';
        pts_in_data = trm_data_all{1}.patients;
    elseif c == 2 % depths
        data = trm_data_all{2}.depths;
        dists = depthsdist;
        comp_name = 'Depth Electrodes';
        pts_in_data = trm_data_all{2}.patients;
    else % strips
        data = trm_data_all{3}.strips;
        dists = stripsdist;
        comp_name = 'Strip Electrodes';
        pts_in_data = trm_data_all{3}.patients;
    end
    
    all_patient_means = cell(length(dists), 1);
    
    % For each distance, collect mean across channels for each patient
    for j = 1:length(dists)
        patient_means = [];

        for p = 1:length(valid_pts)
            pt_idx = find(strcmp(pts_in_data, valid_pts{p}));
            if ~isempty(pt_idx) && j <= size(data,1) && ~isempty(data{j,pt_idx}) && ...
                    size(data{j,pt_idx}, 1) > 0 && size(data{j,pt_idx}, 2) > 0

                patient_mean = nanmean(data{j,pt_idx}, 1); 
                patient_means = [patient_means; patient_mean];
            end
        end

        all_patient_means{j} = patient_means;
        
        if ~isempty(patient_means)
            grand_mean = mean(patient_means, 1); 
            

            plot(frx, grand_mean, '-', 'LineWidth', 3.4, 'Color', cm(max([1 dists(j)]),:));
            hold on;
        end
    end
    

    if c == 1 % grids
        legend({'referential', [num2str(gridsdist(2)) ' mm'], [num2str(gridsdist(3)) ' mm'], [num2str(gridsdist(4)) ' mm'], [num2str(gridsdist(5)) ' mm'], [num2str(gridsdist(6)) ' mm']}, 'location','sw');
    elseif c == 2 % depths
        legend({'referential', [num2str(depthsdist(2)) ' mm'], [num2str(depthsdist(3)) ' mm'], [num2str(depthsdist(4)) ' mm'], [num2str(depthsdist(5)) ' mm']}, 'location','sw');
    else % strips
        legend({'referential', [num2str(stripsdist(2)) ' mm'], [num2str(stripsdist(3)) ' mm']}, 'location','sw');
    end
    
    title(comp_name, 'FontWeight', 'bold', 'FontSize', 12);
    xlabel('Frequency (Hz)');
    ylabel('sqrt(power), rebased');
    set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl);
    grid on;
    set(gca, 'GridAlpha', 0.7, 'MinorGridAlpha', 0.7);
    

    x_50_norm = (log10(50) - log10(2)) / (log10(200) - log10(2));
    x_200_norm = (log10(200) - log10(2)) / (log10(200) - log10(2));
    

    inset_left = pos(1) + pos(3) * x_50_norm;
    inset_width = pos(3) * (x_200_norm - x_50_norm);
    inset_height = pos(4) * 0.70; 
    inset_bottom = pos(2) + pos(4) * 0.25;  
    
    ax_inset = axes('Position', [inset_left, inset_bottom, inset_width, inset_height]);
    

    for j = 1:length(dists)
        patient_means = all_patient_means{j};
        if ~isempty(patient_means)
            grand_mean = mean(patient_means, 1); 
            
            plot(frx, grand_mean, '-', 'LineWidth', 2, 'Color', cm(max([1 dists(j)]),:));
            hold on;
        end
    end
    
    set(gca,'xlim',[50 200],'xscale','log','xtick',[50 100 200],'FontSize',7);
    grid on;
    box on;
end

sgtitle('Aggregated: Mean Across All Patients', 'FontWeight', 'bold', 'FontSize', 14);












