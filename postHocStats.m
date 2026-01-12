% Updated script to do statistics post hoc data analysis
%% Load in data

%function postHocStats

savePostHocPlots = 1;
run2x2plots = 1;
saveSpikeStats = 0;
savePlots = 0;

data_root = getenv("BIPOLAR_DATA");
data_results = getenv("RESULTS");
folderDataBase = data_results;
saveName = {'LL20','LL40','LL100','absDer'};
statsStruct = struct;

numch_eachcond = [];
numchsub_eachcond = [];

outputTable = table;

for jj = 1:length(saveName) 
    processedInt = saveName{jj};
    %folderFigures = folderFiguresCell{jj};
    dataFile  = fullfile(folderDataBase,processedInt,[processedInt '.mat']);
    folderFigures = fullfile(folderDataBase,processedInt);
    load(dataFile);
    
    
    subjCell = [];
    for jjj = 1:length(permResultsCell)
        subjCell{end+1} = permResultsCell{jjj}.subj;
    end
    
    numEls = length(permResultsCell);
    c = brewermap(numEls,'Set3');
    fprintf('numEls: %d\n', numEls);
   
    numChannelsVec = [];
    meanWidthVec = [];
    meanWidthSID = {};
    numChannelsVecSub = [];
    meanWidthVecSub = [];
    numChannelsSID = {};
    subjSpecific = {};
    
    for jjj = 1:length(permResultsCell)
        
        width = permResultsCell{jjj}.meanWidth;
        widthSub = permResultsCell{jjj}.meanWidthSub;
        num = permResultsCell{jjj}.numSigChannels;
        numSub = permResultsCell{jjj}.numSigChannelsSub;
        subj = permResultsCell{jjj}.subj;
        
        if ~isnan(width) & ~isnan(widthSub)
            meanWidthVec = [meanWidthVec width];
            meanWidthVecSub = [meanWidthVecSub widthSub];
            meanWidthSID{end+1} = subj;
            subjSpecific{end+1} = subj;
        end
        
        if ~isnan(num) & ~isnan(numSub)
            numChannelsVec = [numChannelsVec num];
            numChannelsVecSub = [numChannelsVecSub numSub];
            numChannelsSID{end+1} = subj;
        end
        
        
    end

    %% PLOT
    
    if run2x2plots

        colorHist = cmocean('thermal',2);
        colorParallel = brewermap(numEls,'Set3');
        
        histFig = figure;
        tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
        histFig.Position = [839 109 1408 1229];
        nexttile
        histogram(meanWidthVec,BinWidth=10)
        hold on
        histogram(meanWidthVecSub,BinWidth=10)
        legend({'High density','Subsampled'})
        title([processedInt ' Mean Width Histogram'])
        xlim([0 max([meanWidthVec(:);meanWidthVecSub(:)])+10])
        set(gca,'FontSize',14)
        
        nexttile
        histogram(numChannelsVec,BinWidth=1)
        hold on
        histogram(numChannelsVecSub,BinWidth=1)
        legend({'High density','Subsampled'})
        title([processedInt ' Number of Channels Histogram'])
        set(gca,'FontSize',14)
        xlim([0 max([numChannelsVec(:);numChannelsVecSub(:)])+3])

        %How many windows

        numch_eachcond = [numch_eachcond sum(numChannelsVec)];
        numchsub_eachcond = [numchsub_eachcond sum(numChannelsVecSub)];
    
        colormapSpecificWidth = [];
        for jjj = 1:length(subjSpecific)
            ind = find(strcmp(subjCell,subjSpecific(jjj)));
            colormapSpecificWidth = [colormapSpecificWidth; colorParallel(ind,:)];
        end
        
        fprintf('Length HighDensity MeanWidths: %d, Length SubSampled MeanWidths: %d\n', length(meanWidthVec), length(meanWidthVecSub));
        
        ax= nexttile;
        grid(ax,'on')
        xlim([0.5 2.5])
        ylim([0 max([meanWidthVec(:);meanWidthVecSub(:)])+10])
        hold on

        combinedData = [meanWidthVec; meanWidthVecSub];
        group = [ones(size(meanWidthVec)); 2*ones(size(meanWidthVecSub))];
        data1 = meanWidthVec; 
        data2 = meanWidthVecSub;
        vs = violinplot(combinedData, group, 'ShowData', false);
        set(gca, 'XTick', [1 2], 'XTickLabel', {'High Density', 'Subsampled'});
        hold on;
        for i = 1:length(data1)
            plot([1 2], [data1(i) data2(i)], 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
        end
        scatter(ones(size(data1)), data1, 'filled', 'MarkerFaceColor', 'b');
        scatter(2*ones(size(data2)), data2, 'filled', 'MarkerFaceColor', [0.85 0.32 0.1]);
        hold off;

        xlabel('High Density vs. Subsampled')
        ylabel('Mean width')
        title([processedInt ' Mean Spike Width'])
        set(gca,'FontSize',14)
        set(gca, 'XGrid', 'off')
        ax.GridAlpha = 0.4;
        ax.GridColor = [0 0 0];
        
        fprintf('Length HighDensity Channels: %d, Length SubSampled Channels: %d\n', length(numChannelsVec), length(numChannelsVecSub));
        
        ax = nexttile;
        xlim([0.5 2.5])
        grid(ax,'on')
        ylim([0 max([log(numChannelsVec(:)); log(numChannelsVecSub(:))])+1])
        hold on
        set(gca,'FontSize',14)
        %set(gca, 'YTicks', [])
       
        combinedData = [log(numChannelsVec+1); log(numChannelsVecSub+1)];
        group = [ones(size(numChannelsVec)); 2*ones(size(numChannelsVecSub))];
        data1 = log(numChannelsVec+1); 
        data2 = log(numChannelsVecSub+1);
        vs = violinplot(combinedData, group, 'ShowData', false);
        hold on;
        for i = 1:length(data1)
            plot([1 2], [data1(i) data2(i)], 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
        end
        scatter(ones(size(data1)), data1, 'filled', 'MarkerFaceColor', 'b');
        scatter(2*ones(size(data2)), data2, 'filled', 'MarkerFaceColor', [0.85 0.32 0.1]);
        hold off;

        %set(gca, 'YTicks', [])
        title([processedInt ' Mean Spike Width'])
        set(gca,'FontSize',10)
        ax.GridAlpha = 0.4;
        ax.GridColor = [0 0 0];
      
        ylabel('ln(# of Channels)')
        title([processedInt ' Spike Channels'])
        set(gca,'FontSize',10)
        set(gca, 'YTick', [0 1 2 3 4 5])
        set(gca, 'XGrid', 'off')
        set(gca, 'XTick', [1 2], 'XTickLabel', {'', 'High Density', 'Sub-sampled', ''});
        hold on;
        ax.XTickLabel ={'','High density','','Sub-sampled',''};
        ax.GridAlpha = 0.4;
        ax.GridColor = [0 0 0];
        
    
        if savePostHocPlots
            exportgraphics(histFig,fullfile(folderFigures,[processedInt '_2x2.png']),'Resolution',600)
            exportgraphics(histFig,fullfile(folderFigures,[processedInt '_2x2.eps']))
        end
    end

    %%
    histFig = figure;
    tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
    histFig.Position = [800 205 1408 653];
    ax= nexttile;
    grid(ax,'on')
    xlim([0.5 2.5])
    ylim([0 max([meanWidthVec(:);meanWidthVecSub(:)])+10])
    hold on

    combinedData = [meanWidthVec; meanWidthVecSub];
    group = [ones(size(meanWidthVec)); 2*ones(size(meanWidthVecSub))];
    data1 = meanWidthVec;
    data2 = meanWidthVecSub;
    vs = violinplot(combinedData, group, 'ShowData', false);
    set(gca, 'XTick', [1 2], 'XTickLabel', {'High Density', 'Subsampled'});
    hold on;
    for i = 1:length(data1)
        plot([1 2], [data1(i) data2(i)], 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    end
    scatter(ones(size(data1)), data1, 'filled', 'MarkerFaceColor', 'b');
    scatter(2*ones(size(data2)), data2, 'filled', 'MarkerFaceColor', [0.85 0.32 0.1]);
    hold off;

    ylabel('Mean width')
    title([' Mean Spike Width'])
    set(gca,'FontSize',18)
    set(gca, 'XGrid', 'off')
    ax.GridAlpha = 0.4;
    ax.GridColor = [0 0 0];

    fprintf('Length HighDensity Channels: %d, Length SubSampled Channels: %d\n', length(numChannelsVec), length(numChannelsVecSub));

    ax = nexttile;
    xlim([0.5 2.5])
    grid(ax,'on')
    ylim([0 max([log(numChannelsVec(:)); log(numChannelsVecSub(:))])+1])
    hold on

    combinedData = [log(numChannelsVec+1); log(numChannelsVecSub+1)];
    group = [ones(size(numChannelsVec)); 2*ones(size(numChannelsVecSub))];
    data1 = log(numChannelsVec+1);
    data2 = log(numChannelsVecSub+1);
    vs = violinplot(combinedData, group, 'ShowData', false);
    hold on;
    for i = 1:length(data1)
        plot([1 2], [data1(i) data2(i)], 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    end
    scatter(ones(size(data1)), data1, 'filled', 'MarkerFaceColor', 'b');
    scatter(2*ones(size(data2)), data2, 'filled', 'MarkerFaceColor', [0.85 0.32 0.1]);
    hold off;

    title([' Mean Spike Width'])
    ax.GridAlpha = 0.4;
    ax.GridColor = [0 0 0];
    ax.XTick = [0];
    ylabel('natural log (# of Channels)')
    title([' Spike Channels'])
    set(gca, 'YTick', [0 1 2 3 4 5])
    set(gca, 'XGrid', 'off')
    set(gca, 'XTick', [1 2], 'XTickLabel', {'', 'High Density', 'Sub-sampled', ''});
    hold on;
    ax.XTickLabel ={'High density','Sub-sampled','',''};
    ax.GridAlpha = 0.4;
    ax.GridColor = [0 0 0];
    set(gca,'FontSize',18)

    if savePostHocPlots
        exportgraphics(histFig,fullfile(folderFigures,[processedInt '_violon.png']),'Resolution',600)
        exportgraphics(histFig,fullfile(folderFigures,[processedInt '_violin.eps']))
    end

%% COMPUTE STATISTICS

    [pNum,hNum,statsNum] = signrank(numChannelsVec,numChannelsVecSub);
    
    [pWidth,hWidth,statsWidth] = signrank(meanWidthVec,meanWidthVecSub);
    
    statsStruct.name{jj} = processedInt;
    statsStruct.meanWidthSID{jj} = meanWidthSID;
    statsStruct.numChannelsSID{jj} = numChannelsSID;
    statsStruct.pNum{jj} = pNum;
    statsStruct.pWidth{jj} = pWidth;
    statsStruct.meanWidthVec{jj} = meanWidthVec;
    statsStruct.meanWidthVecSub{jj} = meanWidthVecSub;
    statsStruct.numChannelsVec{jj} = numChannelsVec;
    statsStruct.numChannelsVecSubj{jj} = numChannelsVecSub;
    statsStruct.meanWidthDiff{jj} = meanWidthVec - meanWidthVecSub;
    statsStruct.numChannelsDiff{jj} = numChannelsVec - numChannelsVecSub;

end

%% SWARM PLOTS

figure
tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
grid on
hold on
for jj = 1:length(saveName)
    subjSpecific = statsStruct.meanWidthSID{jj};
    colormapSpecific = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecific = [colormapSpecific; c(ind,:)];
    end
    swarmchart(jj,statsStruct.meanWidthDiff{jj},[],colormapSpecific,'filled');
end
title('Differences in Mean Spike Duration by Condition')
xticks([0 1 2 3 4 5])
xticklabels({'','LL20','LL40','LL100','Absolute Derivative',''})

nexttile
grid on
hold on
for jj = 1:length(saveName)
    subjSpecific = statsStruct.numChannelsSID{jj};
    colormapSpecific = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecific = [colormapSpecific; c(ind,:)];
    end
    swarmchart(jj,statsStruct.numChannelsDiff{jj},[],colormapSpecific,'filled');
end
title('Differences in Number of Channels with Spikes by Condition')
xticks([0 1 2 3 4 5])
xticklabels({'','LL20','LL40','LL100','Absolute Derivative',''})
xlabel('condition')

tempFig = gcf;
if savePlots
    exportgraphics(tempFig,fullfile(folderDataBase,['swarm_across_conditions.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderDataBase,['swarm_across_conditions.eps']))
end

%% ACROSS CONDITION PLOTS

allCondData = [];
allCondData_mw = [];
statsCell = {};
figure
tiledlayout(2,1,'TileSpacing','Compact');
ax = nexttile;
grid(ax,'on')

hold on

for jj = 1:length(saveName)
    subjSpecific = statsStruct.meanWidthSID{jj};
  
    colormapSpecific = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecific = [colormapSpecific; c(ind,:)];
        statsCell{ind}.cmap = c(ind,:);
        if jj == 1 
            statsCell{ind}.data{1} = statsStruct.meanWidthDiff{jj}(jjj);
            statsCell{ind}.cond{1} = jj;
        else
            statsCell{ind}.data{end+1} = statsStruct.meanWidthDiff{jj}(jjj);
            statsCell{ind}.cond{end+1} = jj;
                        
        end
    end
end

for jjj =  1:length(subjCell)
    if ~isempty(statsCell{jjj})
        temp = [];
        linePlotPostHoc(jjj) = plot(cell2mat(statsCell{jjj}.cond),cell2mat(statsCell{jjj}.data),'-o','linewidth',4);
    end
end
colororder(colormapSpecific);
for colorInd = 1:length(subjCell)
    if ~isempty(statsCell{colorInd})
        linePlotPostHoc(colorInd).MarkerFaceColor = c(colorInd,:);
    end
end

title('Differences in Mean Spike Duration by Condition')
ylabel('Difference in Mean Spike Duration (ms)')
xticks([0 1 2 3 4 5])
xlim([0 5])
xticklabels({'','','','','',''})
ax = nexttile;
grid(ax,'on')
hold on
for jj = 1:length(saveName)
    subjSpecific = statsStruct.numChannelsSID{jj};
  
    colormapSpecificN = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecificN = [colormapSpecificN; c(ind,:)];
        statsCell{ind}.cmapN = c(ind,:);
        if jj == 1 
            statsCell{ind}.dataN{1} = statsStruct.numChannelsDiff{jj}(jjj);
            statsCell{ind}.condN{1} = jj;
        else
            statsCell{ind}.dataN{end+1} = statsStruct.numChannelsDiff{jj}(jjj);
            statsCell{ind}.condN{end+1} = jj;
                        
        end
    end
end

for jjj =  1:length(subjCell)
    if ~isempty(statsCell{jjj})
        linePlotPostHocN(jjj) = plot(cell2mat(statsCell{jjj}.condN),cell2mat(statsCell{jjj}.dataN),'-o','linewidth',4);
       
    end
end
colororder(colormapSpecificN);
for colorInd = 1:length(subjCell)
    if ~isempty(statsCell{colorInd})
        linePlotPostHocN(colorInd).MarkerFaceColor = c(colorInd,:);
    end
end

title('Differences in Number of Spikes Detected by Condition')
ylabel('Number of channels')
xlabel('Condition')
xticks([0 1 2 3 4 5])
xlim([0 5])
xticklabels({'','LL20','LL40','LL100','Absolute Derivative',''})

tempFig = gcf;
tempFig.Position = [986 636 581 702];

if savePlots
    exportgraphics(tempFig,fullfile(folderDataBase,['across_conditions.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderDataBase,['across_conditions.eps']))
end

%% ACROSS CONDITION PLOTS v2

allCondData = [];
allCondData_mw = [];
statsCell = {};
figure
tiledlayout(2,1,'TileSpacing','Compact');
ax = nexttile;
grid(ax,'on')

hold on

for jj = 1:length(saveName)
    subjSpecific = statsStruct.meanWidthSID{jj};
  
    colormapSpecific = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecific = [colormapSpecific; c(ind,:)];
        statsCell{ind}.cmap = c(ind,:);
        if jj == 1 
            statsCell{ind}.data{1} = statsStruct.meanWidthDiff{jj}(jjj);
            statsCell{ind}.cond{1} = jj;
        else
            statsCell{ind}.data{end+1} = statsStruct.meanWidthDiff{jj}(jjj);
            statsCell{ind}.cond{end+1} = jj;
                        
        end
    end
end

for jjj =  1:length(subjCell)
    if ~isempty(statsCell{jjj})
        temp = [];
        linePlotPostHoc(jjj) = plot(cell2mat(statsCell{jjj}.cond),cell2mat(statsCell{jjj}.data),'-o','linewidth',4);
    end
end
colororder(colormapSpecific);
for colorInd = 1:length(subjCell)
    if ~isempty(statsCell{colorInd})
        linePlotPostHoc(colorInd).MarkerFaceColor = c(colorInd,:);
    end
end

title('Differences in IED Duration')
ylabel('mean spike width (ms)')
xticks([0 1 2 3 4 5])
xlim([0 5])
xticklabels({'','','','','',''})
ax = nexttile;
grid(ax,'on')
hold on
for jj = 1:length(saveName)
    subjSpecific = statsStruct.numChannelsSID{jj};
  
    colormapSpecificN = [];
    for jjj = 1:length(subjSpecific)
        ind = find(strcmp(subjCell,subjSpecific(jjj)));
        colormapSpecificN = [colormapSpecificN; c(ind,:)];
        statsCell{ind}.cmapN = c(ind,:);
        if jj == 1 
            statsCell{ind}.dataN{1} = statsStruct.numChannelsDiff{jj}(jjj);
            statsCell{ind}.condN{1} = jj;
        else
            statsCell{ind}.dataN{end+1} = statsStruct.numChannelsDiff{jj}(jjj);
            statsCell{ind}.condN{end+1} = jj;
                        
        end
    end
end

for jjj =  1:length(subjCell)
    if ~isempty(statsCell{jjj})
        linePlotPostHocN(jjj) = plot(cell2mat(statsCell{jjj}.condN),cell2mat(statsCell{jjj}.dataN),'-o','linewidth',4);
       
    end
end
colororder(colormapSpecificN);
for colorInd = 1:length(subjCell)
    if ~isempty(statsCell{colorInd})
        linePlotPostHocN(colorInd).MarkerFaceColor = c(colorInd,:);
    end
end

title('Differences in Number of Spikes Detected by Condition')
ylabel('Number of channels')
xlabel('Condition')
xticks([0 1 2 3 4 5])
xlim([0 5])
xticklabels({'','LL20','LL40','LL100','Absolute Derivative',''})

tempFig = gcf;
tempFig.Position = [986 636 581 702];

if savePlots
    exportgraphics(tempFig,fullfile(folderDataBase,['across_conditions_v2.png']),'Resolution',600)
    exportgraphics(tempFig,fullfile(folderDataBase,['across_conditions_v2.eps']))
end

%% Save Statistics

if saveSpikeStats
    saveNameSpecific = 'spikeStats.mat';
    save(fullfile(folderDataBase,saveNameSpecific),'statsStruct');
end

