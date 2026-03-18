function violinplot_with_lines(data, labels)

    if isnumeric(data)
        numConds = size(data, 2);
        data_cell = cell(1, numConds);
        for g = 1:numConds
            col = data(:, g);
            data_cell{g} = col(~isnan(col));
        end
    elseif iscell(data)
        numConds = numel(data);
        data_cell = data;
    else
        error('violinplot_with_lines: data must be a numeric matrix or cell array');
    end

    hold on;

    nSubj = max(cellfun(@numel, data_cell));

    for s = 1:nSubj
        % Collect this subject's value across conditions (NaN if missing)
        y = NaN(1, numConds);
        for g = 1:numConds
            if numel(data_cell{g}) >= s
                y(g) = data_cell{g}(s);
            end
        end

        if sum(~isnan(y)) >= 2
            plot(1:numConds, y, '-', ...
                'Color', [0.77 0.77 0.77], 'LineWidth', 3);
        end
    end

    for g = 1:numConds
        x = double(data_cell{g});
        x = x(~isnan(x));

        if isempty(x)
            continue
        end

        [f, xi] = ksdensity(x);
        f = f / max(f) * 0.3;

        fill([g + f, fliplr(g - f)], [xi, fliplr(xi)], ...
             [0.3 0.3 0.3], ...
             'EdgeColor', [0.3 0.3 0.3], ...
             'LineWidth', 2, ...
             'FaceAlpha', 0.4);

        plot(g, median(x), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
        plot([g - 0.1, g + 0.1], [mean(x), mean(x)], 'w-', 'LineWidth', 2.5);
    end

    set(gca, 'XTick', 1:numConds, 'XTickLabel', labels, 'FontSize', 13);
    xlim([0.5, numConds + 0.5]);
    box on;
    grid on;
    ax = gca;
    ax.GridAlpha = 0.5;
end

%{
function violinplot_with_lines(data, labels)
% Custom violin plot function with subject lines and consistent styling.
% 
% Inputs:
%   - data:   [nSubjects x nConditions] matrix
%   - labels: cell array of condition labels (1 x nConditions)

    hold on;
    numConds = size(data, 2);
    nSubj = size(data, 1);

    % --- Plot subject lines first (behind violins) ---
    for s = 1:nSubj
        plot(1:numConds, data(s,:), '-', ...
            'Color', [0.77 0.77 0.77], 'LineWidth', 3);  % gray lines, thicker
    end

    % --- Plot violins ---
    for g = 1:numConds
        x = double(data(:, g));       % ensure type is double
        x = x(~isnan(x));             % remove NaNs

        if isempty(x)
            continue                  % skip if no data
        end

        % Kernel density estimate
        [f, xi] = ksdensity(x);
        f = f / max(f) * 0.3;         % scale width

        % Draw violin with border and transparency
        fill([g + f, fliplr(g - f)], [xi, fliplr(xi)], ...
             [0.3 0.3 0.3], ...                 % gray fill
             'EdgeColor', [0.3 0.3 0.3], ... % dark border
             'LineWidth', 2, ...
             'FaceAlpha', 0.4);               % more transparent

        % Median dot and mean bar
        plot(g, median(x), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 6);
        plot([g - 0.1, g + 0.1], [mean(x), mean(x)], 'w-', 'LineWidth', 2.5);
    end

    % --- Final axis formatting ---
    set(gca, 'XTick', 1:numConds, 'XTickLabel', labels, 'FontSize', 13);
    xlim([0.5, numConds + 0.5]);
    box on;
    grid on;
    ax = gca;
    ax.GridAlpha = 0.5;

end
%}
