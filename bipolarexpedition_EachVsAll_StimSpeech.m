
function [mDiff,mb_m,mARb_m,binz,frx]=bipolarexpedition_EachVsAll_StimSpeech(pt,nchtocheck,windowstocheck)
% BIPOLAR PAIR ANALYSIS: EACH VS. ALL
% see loopbipolarexpedition.m to loop across patients and analyze
% EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)

savePlots = true;
pts = {'EC175','EC183'};

% run on EC175 AND EC183

if ~exist('pt','var')||isempty(pt); pt='EC175'; end %pt='EC175'; % EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)
if ~exist('nchtocheck','var')||isempty(nchtocheck); nchtocheck=128*2; end
if ~exist('windowstocheck','var')||isempty(windowstocheck); windowstocheck=100; end %each window is 1 second of data (non-overlapping)

none1sqrt2log3=3; % 1: no transform, 2: square root, 3: log
g1s2d3=1; % use either grids (1) or strips (2) or depths (3) but not the others
binsz=3; % bin size in mm

xldist=[0 70];
doanglerange=0;
sizeoffont=12;

cm=cool(6); cm(1,:)=[0 0 0];
datadir=getenv("BIPOLAR_DATA");
ptdatadir=fullfile(datadir,'baseline-high-density-data/');
u=dir(ptdatadir); uptbl={}; for i=1:length(u); uname=u(i).name; uptbl{i,1}=uname(1:end-28); end; uptbl(1:2)=[]; clear i u uname
folderFigures = fullfile(datadir,'/results'); if ~exist(folderFigures); mkdir(folderFigures); end


load([datadir '/taggedspikes_April2022']);
sfx=512;
frxrange=[2 200]; %frequency range to examine
ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots

p=find(strcmpi(pts,pt)); %patient number ID
pblocks=strfind(uptbl,pts{p});
for i=1:length(pblocks);
    isbl(i,1)=~isempty(pblocks{i});
end
ptbl=find(isbl); if ~isempty(ptbl); disp(['Loading ' pts{p} ' blocks...']); end

%make vector of stimuli/speech
hasStimTotal = [];
hasSpeechTotal = [];

% load all blocks for this patient and stack their baseline windows together
d=[]; nwind=0;
for b=1:length(ptbl); disp(uptbl{ptbl(b)})
    % load using using "_jk" versions of baseline windows, updated 2/2022
    ptpath = fullfile(ptdatadir, [uptbl{ptbl(b)} '_baselineWindows_fromraw.mat']);
    load(ptpath);
    % get rid of baseline windows containing spikes or artifact
    spksarti=hasspk | hasarti;
    nonspks_windows(:,spksarti)=[];
    hasstim(spksarti)=[]; %update indices for which windows overlap with stimuli/speech
    hasspeech(spksarti)=[];
    hasStimTotal = [hasStimTotal hasstim];
    hasSpeechTotal = [hasSpeechTotal hasspeech];

    clear hasspkvec hasspk hasartivec hasarti spksarti % now clear spike- and artifact-related variables from workspace

    % convert to 3D matrix, combine all windows from consecutive blocks for each patient
    for i=1:size(nonspks_windows,2);
        d(:,:,i+nwind)=nonspks_windows{2,i}';
    end
    nwind=size(d,3);

    clear nonspks_windows info
end; clear b

hasStimTotal = logical(hasStimTotal);
hasSpeechTotal = logical(hasSpeechTotal);

nch=size(d,2);

% load electrode component infor (grid/strip/depth and how many linear contacts they have in a row
% [bpN,bpT]=xlsread(['/Volumes/KLEEN_DRIVE/David/Bipolar project/AN_ElectrodeInfoTDT.xlsx'],pts{p});
% [bpN,bpT]=xlsread(['/Volumes/KLEEN_DRIVE/bipolar_expedition/AN_ElectrodeInfoTDT.xlsx'],pts{p});
an_electrode_info_path = fullfile(datadir, 'AN_ElectrodeInfoTDT.xlsx');
[bpN,bpT]=xlsread(an_electrode_info_path, pts{p});

[em,eleclabels,anatomy]=getelecs(pts{p},2);

% %% if wanting to only look at grids, strips, or depths, then nan the others
% if onlygrids||onlystrips||onlydepths;
%     for r=1:size(bpT,1)
%         if any(strcmpi(bpT(r,2),{'grid','minigrid'})) && ~onlygrids  || ...
%                 strcmpi(bpT(r,2),'strip')              && ~onlystrips || ...
%                 strcmpi(bpT(r,2),'depth')              && ~onlydepths;
%             d(:,bpN(r,1):bpN(r,2),:)=nan;
%         end
%     end; clear r
% end

%% look at either grids, strips, or depths, and nan the others
          for r=1:size(bpT,1)
            if [g1s2d3~=1 && any(strcmpi(bpT(r,2),{'grid','minigrid'}))]  || ...
               [g1s2d3~=2 &&     strcmpi(bpT(r,2),'strip')]               || ...
               [g1s2d3~=3 &&     strcmpi(bpT(r,2),'depth')];
                       d(:,bpN(r,1):bpN(r,2),:)=nan; %nan out component that isn't relevant to this run (see g1s2d3 above)
              bp_distance(bpN(r,1):bpN(r,2))  =nan; %nan out their corresponding distances (irrelevant for this run)
              bp_angle   (bpN(r,1):bpN(r,2))  =nan; %and angles, similarly
            end
          end; clear r


%% bad channels
badchI=isnan(mean(mean(d,1),3))'; %any channel with nans in any window will be a "bad channel"
badchI(badchidx{p})=1; %follow up with all marked as bad channels in original preprocessing
x=[]; xbch=true(size(badchI)); for i=1:size(bpT,1); x=[x bpN(i,1):bpN(i,2)]; end; xbch(x)=false;
badchI(xbch)=1; %find and nan empty channels unaccounted for by component rows
d(:,badchI,:)=nan;
okc=~badchI; clear x xbch


%*at this point, isolate speech or stim or non-speech/stim

%*at this point, isolate STG or IFG channels
[STG]=getelecs_region(pt,'stg',2);
[IFG]=getelecs_region(pt,{'po','pt'},2);

windowstocheck=min([windowstocheck size(d,3)]);
windowstocheck=1:windowstocheck; %convert to a vector of windows, 1:X


%% ALL PAIRS (each vs. all others) analysis and example plot
% ***hint hint: opportunity here to select speech or stim windows!

clear Straces_allch;

dSpeech = d(:,:,(hasSpeechTotal & ~hasStimTotal));
dNoSpeechNoStim = d(:,:,(~hasSpeechTotal & ~hasStimTotal));
dStim = d(:,:,(hasStimTotal & ~hasSpeechTotal));

%NEW BELOW 
%dSpeech =           dSpeech         (:,:,1:min([size(d,3) size(dSpeech,3)]));
%dNoSpeechNoStim =   dNoSpeechNoStim (:,:,1:min([size(d,3) size(dNoSpeechNoStim,3)]));
%dStim =             dStim           (:,:,1:min([size(d,3) size(dStim,3)]));

% ORIGINAL BELOW
dSpeech = dSpeech(:,:,windowstocheck);
dNoSpeechNoStim = dNoSpeechNoStim(:,:,windowstocheck);
dStim = dStim(:,:,windowstocheck);


disp('Speech windows'); [mSpeech,~,~,~]=bpspectra_EachVsAll_2025(dSpeech,sfx,frxrange,em,nchtocheck,none1sqrt2log3);
disp('NoSpeechNoStim windows'); [mNoSpeechNoStim,~,~,~]=bpspectra_EachVsAll_2025(dNoSpeechNoStim,sfx,frxrange,em,nchtocheck,none1sqrt2log3);
disp('Stim windows'); [mStim,~,~,~]=bpspectra_EachVsAll_2025(dStim,sfx,frxrange,em,nchtocheck,none1sqrt2log3);

%mSpeech = M(:,:,:,(hasSpeechTotal & ~hasStimTotal));
%mNoSpeechNoStim = M(:,:,:,(~hasSpeechTotal & ~hasStimTotal));
%mStim = M(:,:,:,(hasStimTotal & ~hasSpeechTotal));

d=d(:,:,windowstocheck); clear Straces_allch; %free up RAM by getting rid of whatever won't be used (only using first ___ number of windows)
[M,Mrefave,Mbpdist,frx]=bpspectra_EachVsAll_2025(d,sfx,frxrange,em,nchtocheck,none1sqrt2log3);

%% mean across windows

M=sq(mean(M,4));

%% now that we have bipolar pairs and spectra, we can select the channels we are interested in

% which channnels are we interseted in
chansIntSTG = STG;
chansIntSTG = chansIntSTG(chansIntSTG<=256);

mSpeechSTG = mSpeech(chansIntSTG,chansIntSTG,:,:);
mStimSTG = mStim(chansIntSTG,chansIntSTG,:,:);
mNoSpeechNoStimSTG = mNoSpeechNoStim(chansIntSTG,chansIntSTG,:,:);
MbpdistSTG = Mbpdist(chansIntSTG, chansIntSTG);

chansIntIFG = IFG;
chansIntIFG = chansIntIFG(chansIntIFG<=256);

mSpeechIFG = mSpeech(chansIntIFG,chansIntIFG,:,:);
mStimIFG = mStim(chansIntIFG,chansIntIFG,:,:);
mNoSpeechNoStimIFG = mNoSpeechNoStim(chansIntIFG,chansIntIFG,:,:);
MbpdistIFG = Mbpdist(chansIntIFG,chansIntIFG);


%% unstack and line up into 2D matrix to create channels^2 X frequencies for easier indexing-->binning
Mflat=[];
Mflat_bpdist=[];
for c2=1:nchtocheck;
    Mflat=[Mflat; sq(M(:,c2,:))];
    Mflat_bpdist=[Mflat_bpdist; Mbpdist(:,c2)]; % corresponding distance index
end

binz=0:binsz:85; % can change binsz PRN at this point to test different bin resolutions on the plots below
clear mx my mz
nbinz=length(binz)-1;
mz=[];
%bin by distance and take the mean, creating: frequency X distance (binned)
for i=1:nbinz
    mz(:,i)=nanmean(Mflat(Mflat_bpdist>=binz(i) & Mflat_bpdist<binz(i+1),:),1);
    binindex_min_ltnmax(i,:)=[binz(i) binz(i+1)];
end

% Notes:
% mz is frequency X distance (binned), after having averaged across windows
% binindex_min_ltmax tells you for each column of mz what was the
%         1] minimum (>0) distance, and
%         2] the "less than max" (<) distance
%         that were used to index bipolar pairs for that bin
%         Note: first bin includes zero, which corresponds to
%         the bin containing the referential channels


%% unstack and line up into 3D matrix to create channels^2 X frequencies x trials for easier indexing-->binning
% NOTE: the transform (based on none1sqrt2log3) is performed INSIDE this
[mz_zSpeech,mzSpeech,mflatSpeech,bpdistSpeech ] =               bin_zscore_trial(mSpeech,nchtocheck,Mbpdist,binsz,frx,none1sqrt2log3);
[mz_zStim,mzStim,mflatStim,bpdistStim ] =                       bin_zscore_trial(mStim,nchtocheck,Mbpdist,binsz,frx,none1sqrt2log3);
[mz_zNoST,mzNoST,mflatNoST,bpdistNoST ] =                       bin_zscore_trial(mNoSpeechNoStim,nchtocheck,            Mbpdist,   binsz,frx,none1sqrt2log3);

[mz_zSpeechSTG,mzSpeechSTG,mflatSpeechSTG,bpdistSpeechSTG ] =   bin_zscore_trial(mSpeechSTG,length(chansIntSTG),MbpdistSTG,binsz,frx,none1sqrt2log3);
[mz_zStimSTG,mzStimSTG,mflatStimSTG,bpdistStimSTG ] =           bin_zscore_trial(mStimSTG,length(chansIntSTG),MbpdistSTG,binsz,frx,none1sqrt2log3);
[mz_zNoSTSTG,mzNoSTSTG,mflatNoSTSTG,bpdistNoSTSTG ] =           bin_zscore_trial(mNoSpeechNoStimSTG,length(chansIntSTG),MbpdistSTG,binsz,frx,none1sqrt2log3);

[mz_zSpeechIFG,mzSpeechIFG,mflatSpeechIFG,bpdistSpeechIFG ] =   bin_zscore_trial(mSpeechIFG,length(chansIntIFG),MbpdistIFG,binsz,frx,none1sqrt2log3);
[mz_zStimIFG,mzStimIFG,mflatStimIFG,bpdistStimIFG ] =           bin_zscore_trial(mStimIFG,length(chansIntIFG),MbpdistIFG,binsz,frx,none1sqrt2log3);
[mz_zNoSTIFG,mzNoSTIFG,mflatNoSTIFG,bpdistNoSTIFG ] =           bin_zscore_trial(mNoSpeechNoStimIFG,length(chansIntIFG),MbpdistIFG,binsz,frx,none1sqrt2log3);


%%  permutation testing

%for chan = 1:size(mSpeech,1)
[sizeX,~,sizeY] = size(mz_zSpeech);
[clustersSpeech, pValuesSpeech, tSumsSpeech, permutationDistributionSpeech] = permutest(permute(mz_zSpeech,[1,3,2]),permute(mz_zNoST,[1,3,2]),false,[],[],1);
[clustersStim, pValuesStim, tSumsStim, permutationDistributionStim] = permutest(permute(mz_zStim,[1,3,2]),permute(mz_zNoST,[1,3,2]),false,[],[],1);
[clustersStimSpeech, pValuesStimSpeech, tSumsStimSpeech, permutationDistributionStimSpeech] = permutest(permute(mz_zStim,[1,3,2]),permute(mz_zSpeech,[1,3,2]),false,[],[],1);

%for chan = 1:size(mSpeech,1)
[sizeXSTG,~,sizeYSTG] = size(mz_zSpeechSTG);
[clustersSpeechSTG, pValuesSpeechSTG, tSumsSpeechSTG, permutationDistributionSpeechSTG] = permutest(permute(mz_zSpeechSTG,[1,3,2]),permute(mz_zNoSTSTG,[1,3,2]),false,[],[],1);
[clustersStimSTG, pValuesStimSTG, tSumsStimSTG, permutationDistributionStimSTG] = permutest(permute(mz_zStimSTG,[1,3,2]),permute(mz_zNoSTSTG,[1,3,2]),false,[],[],1);
[clustersStimSpeechSTG, pValuesStimSpeechSTG, tSumsStimSpeechSTG, permutationDistributionStimSpeechSTG] = permutest(permute(mz_zStimSTG,[1,3,2]),permute(mz_zSpeechSTG,[1,3,2]),false,[],[],1);

%for chan = 1:size(mSpeech,1)
[sizeXIFG,~,sizeYIFG] = size(mz_zSpeechIFG);
[clustersSpeechIFG, pValuesSpeechIFG, tSumsSpeechIFG, permutationDistributionSpeechIFG] = permutest(permute(mz_zSpeechIFG,[1,3,2]),permute(mz_zNoSTIFG,[1,3,2]),false,[],[],1);
[clustersStimIFG, pValuesStimIFG, tSumsStimIFG, permutationDistributionStimIFG] = permutest(permute(mz_zStimIFG,[1,3,2]),permute(mz_zNoSTIFG,[1,3,2]),false,[],[],1);
[clustersStimSpeechIFG, pValuesStimSpeechIFG, tSumsStimSpeechIFG, permutationDistributionStimSpeechIFG] = permutest(permute(mz_zStimIFG,[1,3,2]),permute(mz_zSpeechIFG,[1,3,2]),false,[],[],1);


%%

[clusterSigSpeech,pValsSigSpeech,boundarySigSpeech,clustXSpeech,clustYSpeech] = signif_boundary(clustersSpeech,pValuesSpeech,sizeX,sizeY);
[clusterSigStim,pValsSigStim,boundarySigStim,clustXStim,clustYStim] = signif_boundary(clustersStim,pValuesStim,sizeX,sizeY);
[clusterSigStimSpeech,pValsSigStimSpeech,boundarySigStimSpeech,clustXStimSpeech,clustYStimSpeech] = signif_boundary(clustersStimSpeech,pValuesStimSpeech,sizeX,sizeY);

[clusterSigSpeechSTG,pValsSigSpeechSTG,boundarySigSpeechSTG,clustXSpeechSTG,clustYSpeechSTG] = signif_boundary(clustersSpeechSTG,pValuesSpeechSTG,sizeXSTG,sizeYSTG);
[clusterSigStimSTG,pValsSigStimSTG,boundarySigStimSTG,clustXStimSTG,clustYStimSTG] = signif_boundary(clustersStimSTG,pValuesStimSTG,sizeXSTG,sizeYSTG);
[clusterSigStimSpeechSTG,pValsSigStimSpeechSTG,boundarySigStimSpeechSTG,clustXStimSpeechSTG,clustYStimSpeechSTG] = signif_boundary(clustersStimSpeechSTG,pValuesStimSpeechSTG,sizeXSTG,sizeYSTG);

[clusterSigSpeechIFG,pValsSigSpeechIFG,boundarySigSpeechIFG,clustXSpeechIFG,clustYSpeechIFG] = signif_boundary(clustersSpeechIFG,pValuesSpeechIFG,sizeXIFG,sizeYIFG);
[clusterSigStimIFG,pValsSigStimIFG,boundarySigStimIFG,clustXStimIFG,clustYStimIFG] = signif_boundary(clustersStimIFG,pValuesStimIFG,sizeXIFG,sizeYIFG);
[clusterSigStimSpeechIFG,pValsSigStimSpeechIFG,boundarySigStimSpeechIFG,clustXStimSpeechIFG,clustYStimSpeechIFG] = signif_boundary(clustersStimSpeechIFG,pValuesStimSpeechIFG,sizeXIFG,sizeYIFG);

%%
avgNoST = nanmean(mz_zNoST,2);
avgSpeech = nanmean(mz_zSpeech,2);
avgStim = nanmean(mz_zStim,2);
maxZ = max([avgNoST(:);avgSpeech(:);avgStim(:)]);
minZ = min([avgNoST(:);avgSpeech(:);avgStim(:)]);
maxAbs = max(abs(maxZ),abs(minZ));

avgNoSTSTG = nanmean(mz_zNoSTSTG,2);
avgSpeechSTG = nanmean(mz_zSpeechSTG,2);
avgStimSTG = nanmean(mz_zStimSTG,2);
maxZSTG = max([avgNoSTSTG(:);avgSpeechSTG(:);avgStimSTG(:)]);
minZSTG = min([avgNoSTSTG(:);avgSpeechSTG(:);avgStimSTG(:)]);
maxSTGAbs = max(abs(maxZSTG),abs(minZSTG));

%%
avgSpeechBase = squeeze(nanmean(mz_zSpeech,2)) - squeeze(nanmean(mz_zNoST,2));
avgStimBase = squeeze(nanmean(mz_zStim,2)) - squeeze(nanmean(mz_zNoST,2));
avgSpeechStim = squeeze(nanmean(mz_zSpeech,2)) - squeeze(nanmean(mz_zStim,2));
maxZsub = max([avgSpeechBase(:);avgStimBase(:);avgSpeechStim(:)]);
minZsub = min([avgSpeechBase(:);avgStimBase(:);avgSpeechStim(:)]);
maxZsubAbs = max(abs(maxZsub),abs(minZsub));

avgNoSTIFG = nanmean(mz_zNoSTIFG,2);
avgSpeechIFG = nanmean(mz_zSpeechIFG,2);
avgStimIFG = nanmean(mz_zStimIFG,2);
maxZIFG = max([avgNoSTIFG(:);avgSpeechIFG(:);avgStimIFG(:)]);
minZIFG = min([avgNoSTIFG(:);avgSpeechIFG(:);avgStimIFG(:)]);
maxIFGabs = max(abs(minZIFG),abs(maxZIFG));

avgSpeechBaseSTG = squeeze(nanmean(mz_zSpeechSTG,2)) - squeeze(nanmean(mz_zNoSTSTG,2));
avgStimBaseSTG = squeeze(nanmean(mz_zStimSTG,2)) - squeeze(nanmean(mz_zNoSTSTG,2));
avgSpeechStimSTG = squeeze(nanmean(mz_zSpeechSTG,2)) - squeeze(nanmean(mz_zStimSTG,2));
maxZsubSTG = max([avgSpeechBaseSTG(:);avgStimBaseSTG(:);avgSpeechStimSTG(:)]);
minZsubSTG = min([avgSpeechBaseSTG(:);avgStimBaseSTG(:);avgSpeechStimSTG(:)]);
maxSTGAbs=max(abs(maxZsubSTG),abs(minZsubSTG));

%
avgSpeechBaseClust = nan(size(avgSpeechBase));
for jj = 1:length(clusterSigSpeech)
    clusterTemp = clusterSigSpeech{jj};
    avgSpeechBaseClust(clusterTemp) = avgSpeechBase(clusterTemp);
end;

avgStimBaseClust = nan(size(avgStimBase));
for jj = 1:length(clusterSigStim)
    clusterTemp = clusterSigStim{jj};
    avgStimBaseClust(clusterTemp) = avgStimBase(clusterTemp);
end;

avgSpeechStimClust = nan(size(avgSpeechStim));
for jj = 1:length(clusterSigStimSpeech)
    clusterTemp = clusterSigStimSpeech{jj};
    avgSpeechStimClust(clusterTemp) = avgSpeechStim(clusterTemp);
end;

%
avgSpeechBaseClustSTG = nan(size(avgSpeechBaseSTG));
for jj = 1:length(clusterSigSpeechSTG)
    clusterTemp = clusterSigSpeechSTG{jj};
    avgSpeechBaseClustSTG(clusterTemp) = avgSpeechBaseSTG(clusterTemp);
end;

avgStimBaseClustSTG = nan(size(avgStimBaseSTG));
for jj = 1:length(clusterSigStimSTG)
    clusterTemp = clusterSigStimSTG{jj};
    avgStimBaseClustSTG(clusterTemp) = avgStimBaseSTG(clusterTemp);
end;

avgSpeechStimClustSTG = nan(size(avgSpeechStimSTG));
for jj = 1:length(clusterSigStimSpeechSTG)
    clusterTemp = clusterSigStimSpeechSTG{jj};
    avgSpeechStimClustSTG(clusterTemp) = avgSpeechStimSTG(clusterTemp);
end;

avgSpeechBaseIFG = squeeze(nanmean(mz_zSpeechIFG,2)) - squeeze(nanmean(mz_zNoSTIFG,2));
avgStimBaseIFG = squeeze(nanmean(mz_zStimIFG,2)) - squeeze(nanmean(mz_zNoSTIFG,2));
avgSpeechStimIFG = squeeze(nanmean(mz_zSpeechIFG,2)) - squeeze(nanmean(mz_zStimIFG,2));
maxZsubIFG = max([avgSpeechBaseIFG(:);avgStimBaseIFG(:);avgSpeechStimIFG(:)]);
minZsubIFG = min([avgSpeechBaseIFG(:);avgStimBaseIFG(:);avgSpeechStimIFG(:)]);
maxIFGabs = max(abs(minZsubIFG),abs(maxZsubIFG));
%
avgSpeechBaseClustIFG = nan(size(avgSpeechBaseIFG));
for jj = 1:length(clusterSigSpeechIFG)
    clusterTemp = clusterSigSpeechIFG{jj};
    avgSpeechBaseClustIFG(clusterTemp) = avgSpeechBaseIFG(clusterTemp);
end;

avgStimBaseClustIFG = nan(size(avgStimBaseIFG));
for jj = 1:length(clusterSigStimIFG)
    clusterTemp = clusterSigStimIFG{jj};
    avgStimBaseClustIFG(clusterTemp) = avgStimBaseIFG(clusterTemp);
end;

avgSpeechStimClustIFG = nan(size(avgSpeechStimIFG));
for jj = 1:length(clusterSigStimSpeechIFG)
    clusterTemp = clusterSigStimSpeechIFG{jj};
    avgSpeechStimClustIFG(clusterTemp) = avgSpeechStimIFG(clusterTemp);
end;
%


%% plot condensed power


mz_zStimSTG_gamma = permute(squeeze(mean((mz_zStimSTG(frx>=50,:,:)),1)),[2,1]);
mz_zNoSTSTG_gamma = permute(squeeze(mean((mz_zNoSTSTG(frx>=50,:,:)),1)),[2,1]);
mz_zStim_gamma = permute(squeeze(mean((mz_zStim(frx>=50,:,:)),1)),[2,1]);
mz_zNoST_gamma = permute(squeeze(mean((mz_zNoST(frx>=50,:,:)),1)),[2,1]);

binzPlotSTG = binz(2:size(mz_zStimSTG_gamma,1)+1);
binzPlotTotal = binz(2:size(mz_zStim_gamma,1)+1);


[clustersSpeechSTG_gamma, pValuesSpeechSTG_gamma, tSumsSpeechSTG_gamma, permutationDistributionSpeechSTG_gamma] = permutest(mz_zStimSTG_gamma,mz_zNoSTSTG_gamma,false,[],[],1);
[clustersSpeech_gamma, pValuesSpeech_gamma, tSumsSpeech_gamma, permutationDistributionSpeech_gamma] = permutest(mz_zStim_gamma,mz_zNoST_gamma,false,[],[],1);


%can save files for plotting later in the fig5_out file
save(fullfile(folderFigures,['/stg_Devon_' pt(3:end) '.mat']),'avgStimBaseSTG','maxSTGAbs', ...
    'dSpeech','mz_zSpeechSTG','mz_zNoSTSTG','dNoSpeechNoStim', 'dStim','mSpeech', 'mNoSpeechNoStim',...
    'avgStimBaseClustSTG','mz_zStimSTG_gamma','mz_zNoSTSTG_gamma','clustersSpeechSTG_gamma','pValuesSpeechSTG_gamma', ...
    'binzPlotSTG','MbpdistSTG','binzPlotTotal','binz','ft','ftl',...
    'mStim', 'Mbpdist', 'frx', '-v7.3');

%% 