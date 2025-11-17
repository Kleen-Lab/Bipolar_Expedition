
% (1) Runs signed rank test on # chn w/ significant spikes (HD vs subsamp)
% (2) Runs signed rank test on mean width (HD vs subsamp)
% (3) Runs ANOVA on difference in # channels w/ spikes by condition
% (4) Runs ANOVA on difference in mean spike duration (ms) by condition

folderDataBase = '/scratch/bipolar_expedition';

spikeStats = load(fullfile(folderDataBase,'spikeStatsV3.mat'));
conds = {'LL20', 'LL40', 'LL100', 'absDer'};


%% (1)

% run the SRT on the pair of HD vs subsampled for each condition

A = [];
Asub = [];
channelsCell = spikeStats.statsStruct.numChannelsVec;
channelsCellSub = spikeStats.statsStruct.numChannelsVecSubj;

for i =1:length(conds)
    A = [A; channelsCell{i}];
    Asub = [Asub; channelsCellSub{i}];
end

A = A';
Asub = Asub';

for i = 1:length(conds)
    [p, h, stats] = signrank(A(:,i), Asub(:,i));
    disp(p);
end


%% (2)

% run the SRT on the pair of HD vs subsampled for each condition


B = [];
Bsub = [];
widthCell = spikeStats.statsStruct.meanWidthVec;
widthCellSub = spikeStats.statsStruct.meanWidthVecSub;

for i=1:length(conds)
%    if i==3 %LL 100
%        vecB = [NaN, widthCell{i}]; 
%        vecBsub = [NaN, widthCellSub{i}];
%        B = [B; vecB];
%        Bsub = [Bsub; vecBsub];
 %       continue;
  %  end
    B = [B; widthCell{i}];
    Bsub = [Bsub; widthCellSub{i}];
end

B = B';
Bsub = Bsub';

% need to remove patient 1 row (nan exists)
B(any(isnan(B), 2), :) = [];
Bsub(any(isnan(Bsub), 2), :) = [];

for i = 1:length(conds)
    [p, h, stats] = signrank(B(:,i), Bsub(:,i));
    disp(p);
end

%% (3)

Achannels = [];
diffChannelsCell = spikeStats.statsStruct.numChannelsDiff;

for i = 1:length(conds)
    Achannels = [Achannels; diffChannelsCell{i}];
end

Achannels = Achannels';
T1 = array2table(Achannels, 'VariableNames', conds);

rm = fitrm(T1, 'LL20-absDer~1', 'WithinDesign', conds);
ranovaTable = ranova(rm);
disp(ranovaTable.pValue);

%% (4)

Awidth = [];
diffWidthsCell = spikeStats.statsStruct.meanWidthDiff;

for i = 1:length(conds)
%    if i==3
%        vec = [NaN, diffWidthsCell{i}];
%        Awidth = [Awidth; vec];
 %       continue;
 %   end
    Awidth = [Awidth; diffWidthsCell{i}];
end

Awidth = Awidth';

% need to remove patient 1 row (nan exists)
Awidth(any(isnan(Awidth), 2), :) = [];

T2 = array2table(Awidth, 'VariableNames', conds);

rm2 = fitrm(T2, 'LL20-absDer~1', 'WithinDesign', conds);
ranovaTable2 = ranova(rm2);
disp(ranovaTable2.pValue);




