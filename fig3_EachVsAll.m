
function [pt, binz, toplot, frx, binsz, Mbp_distance, cm_distance]=fig3_EachVsAll(pt,nchtocheck,windowstocheck)

% script to save data for plotting figure 3 omni-directional BPRR analysis
% for either 'EC175' or 'EC133' (both have intact 16 x 16 grids)
% Run this script to save ECXXX_fig3_data.mat
% Then plot the figure using fig3_plot.m with that .mat file

    data_root = getenv("BIPOLAR_DATA");
    data_results = getenv("RESULTS");

    save_mat = true; % save .mat file for later plotting with fig3_plot.m
    
    if ~exist('pt','var')||isempty(pt); pt='EC175'; end % EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)    
    if ~exist('nchtocheck','var')||isempty(nchtocheck); nchtocheck=128*2; end
    if ~exist('windowstocheck','var')||isempty(windowstocheck); windowstocheck=250; end %each window is 1 second of data (non-overlapping)
    
    none1sqrt2log3=2; % 1: no transform, 2: square root, 3: log
    g1s2d3=1; % use either grids (1) or strips (2) or depths (3) but not the others
    binsz=2; % bin size in mm
    xldist=[0 60];
    doanglerange=0;
    onlygrids=true;
    onlydepths=false;
    onlystrips=false;
    sizeoffont=12;
    
    cm=cool(6); cm(1,:)=[0 0 0];
    
    % baseline-high-density-data folder is and where to save the .mat file to
    datadir = fullfile(data_root, 'baseline-high-density-data');
    %cd([datadir])
    tag_spikes_path = fullfile(data_root, 'taggedspikes_April2022.mat');
    load(tag_spikes_path);
    sfx=512;
    frxrange=[2 200]; %frequency range to examine
      ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots
     
    u=dir(datadir); uptbl={}; for i=1:length(u); uname=u(i).name; uptbl{i,1}=uname(1:end-28); end; uptbl(1:2)=[]; clear i u uname
    
    p=find(strcmpi(pts,pt)); %patient number ID
    pblocks=strfind(uptbl,pts{p});
    for i=1:length(pblocks);
        isbl(i,1)=~isempty(pblocks{i});
    end
    ptbl=find(isbl); if ~isempty(ptbl); disp(['Loading ' pts{p} ' blocks...']); end
    
    % load all blocks for this patient and stack their baseline windows together
    d=[]; nwind=0;
    for b=1:length(ptbl); disp(uptbl{ptbl(b)})
        % load using using "_jk" versions of baseline windows, updated 2/2022
        ptpath = fullfile(datadir, [uptbl{ptbl(b)} '_baselineWindows_fromraw.mat']);
        load(ptpath);
        % get rid of baseline windows containing spikes or artifact
        spksarti=hasspk | hasarti;
        nonspks_windows(:,spksarti)=[];
        hasstim(spksarti)=[]; %update indices for which windows overlap with stimuli/speech
        hasspeech(spksarti)=[];
        clear hasspkvec hasspk hasartivec hasarti spksarti % now clear spike- and artifact-related variables from workspace
    
        % convert to 3D matrix, combine all windows from consecutive blocks for each patient
        for i=1:size(nonspks_windows,2)
            d(:,:,i+nwind)=nonspks_windows{2,i}';
        end
        nwind=size(d,3);
    
        clear nonspks_windows info
    end; clear b
    
    nch=size(d,2);
    
    % load electrode component infor (grid/strip/depth and how many linear contacts they have in a row

    an_electrode_info_path = fullfile(data_root, 'AN_ElectrodeInfoTDT.xlsx');
    [bpN,bpT]=xlsread(an_electrode_info_path, pts{p});
    
    
    [em,~,~]=getelecs(pts{p},2); % call to elec coordinates fxn
    
    cm=cool(6); cm(1,:)=[0 0 0];
 

    sfx=512;
    frxrange=[2 200]; %frequency range to examine
      ft=[2 5 10 20 50 100 200]; ftl=cellstr(num2str(ft')); %frequency labels for plots
     
    
    
    %% look at either grids, strips, or depths, and nan the others
    % here we are only looking at grids
              for r=1:size(bpT,1)
                if [g1s2d3~=1 && any(strcmpi(bpT(r,2),{'grid','minigrid'}))]  || ...
                   [g1s2d3~=2 &&     strcmpi(bpT(r,2),'strip')]               || ...
                   [g1s2d3~=3 &&     strcmpi(bpT(r,2),'depth')];
                           d(:,bpN(r,1):bpN(r,2),:)=nan; %nan out component that isn't relevant to this run (see g1s2d3 above)
                  bp_distance(bpN(r,1):bpN(r,2))  =nan; %nan out their corresponding distances (irrelevant for this run)
                  bp_angle   (bpN(r,1):bpN(r,2))  =nan; %and angles, similarly
                end
              end; clear r
    
    
    % bad channels
    badchI=isnan(mean(mean(d,1),3))'; %any channel with nans in any window will be a "bad channel"
    badchI(badchidx{p})=1; %follow up with all marked as bad channels in original preprocessing
    x=[]; xbch=true(size(badchI)); for i=1:size(bpT,1); x=[x bpN(i,1):bpN(i,2)]; end; xbch(x)=false;
    badchI(xbch)=1; %find and nan empty channels unaccounted for by component rows
    d(:,badchI,:)=nan;
    okc=~badchI; clear x xbch
    
    
    
    windowstocheck=min([windowstocheck size(d,3)]);
    windowstocheck=1:windowstocheck; %convert to a vector of windows, 1:X
    
    %% ALL PAIRS (each vs. all others) analysis and example plot
     d=d(:,:,windowstocheck); clear Straces_allch; %free up RAM by getting rid of whatever won't be used (only using first ___ number of windows)
                     % ***opportunity here to select speech or stim windows
     [M,M_averef,Mbp_distance,frx,~,Mbp_angle]=bpspectra_EachVsAll_2025(d,sfx,frxrange,em,nchtocheck, none1sqrt2log3);
         % M is bipolarchannel1 X bipolarchannel2 X frx X 1secwindow
     nfrx=length(frx);
     
    
    %% Unstack and line up into 2D matrix to create channels^2 X frequencies for easier indexing-->binning
    mb=[];
    mARb=[];
    binz=[-1 0:binsz:85]; % can change binsz PRN at this point to test different bin resolutions on the plots below
    clear mx my mz  
    nbinz=length(binz)-1;
     parfor w=windowstocheck; disp(num2str(w)); % parfor here to run the windows
        Mflat=[];
          Mbp_bpflat_distance=[];
          M_bpflat_angle=[];
        M_averef_flat=[];
        for c2=1:nchtocheck;
            Mflat=[Mflat; squeeze(M(:,c2,:,w))];              
            Mbp_bpflat_distance=[Mbp_bpflat_distance; Mbp_distance(:,c2)]; % corresponding distance index
            M_bpflat_angle   =[M_bpflat_angle;    Mbp_angle(:,c2)]; % corresponding angle index
            M_averef_flat=[M_averef_flat; squeeze(M_averef(:,c2,:,w))];              
        end
       
        %bin by distance and take the mean, creating: frequency X binned distance
           % first for referential (distance = 0)
        for i=1:nbinz % binz will be including >lower bound and up to and including (<=) upper bound
            mb  (:,i,w)=nanmean(Mflat        (Mbp_bpflat_distance>binz(i) & Mbp_bpflat_distance<=binz(i+1),:),1);
            mARb(:,i,w)=nanmean(M_averef_flat(Mbp_bpflat_distance>binz(i) & Mbp_bpflat_distance<=binz(i+1),:),1); %AR means Average Reference
        end
       
        disp([num2str(round(windowstocheck(w)/windowstocheck(end)*100,1)) '% of windows'])
    end; disp('Done')
    
    binz(1)=[]; %remove negative 1, only there to get referential channels (distance=0) for the parfor loop
    
    %         mb(:,1,w)  =nanmean(Mflat      (Mflatbp_distance==0,:),1);
    %         mARb(:,1,w)=nanmean(Mrefaveflat(Mflatbp_distance==0,:),1);
    
    
    %% --- Figure 3 --- (Omnidirectional re-referencing)
    
    % % NOTES:
    % % mb is frequency X binned distance X windows, and you can average across windows
    % % binindex_min_ltmax tells you for each column of mb what was the
    % %         1] minimum (>0) distance, and
    % %         2] the "less than max" (<) distance
    % %         that were used to index bipolar pairs for that bin
    % %         Note: first bin includes zero, which corresponds to
    % %         the bin containing the referential channels
    
    
    % % zscore the log transformed power according to frequency
    
    % transform power before calculating
    switch none1sqrt2log3
      case 2
            mb__m=squeeze(mean(sqrt(mb),3));
            mARb__m=squeeze(mean(sqrt(mARb),3));
             nfrx=length(frx);
             mb__m_z=mb__m; mARb__m_z=mARb__m;
             for i=1:nfrx; nns=~isnan(mb__m(i,:));
                 mb__m_z(i,nns)=zscore((mb__m(i,nns)));
                 mARb__m_z(i,nns)=zscore((mARb__m(i,nns)));
             end;
      case 3
            mblogm=squeeze(mean(log(mb),3));
            mARblogm=squeeze(mean(log(mARb),3));
             nfrx=length(frx);
             mblogm_z=mblogm; mARblogm_z=mARblogm;
             for i=1:nfrx; nns=~isnan(mblogm(i,:));
                 mblogm_z(i,nns)=zscore((mblogm(i,nns)));
                 mARblogm_z(i,nns)=zscore((mARblogm(i,nns)));
             end;
    end
    
    
     figure('color','w','position',[[54 223 1907 1102]]); colormap(parula);
     
     % plot log-transformed version
     toplot=mean((mb__m_z),3);  % mb or mb_z
     cm_distance=flipud(cmocean('deep',1+ceil(max(max(Mbp_distance)))));
    
    % plot confusion matrix of bipolar distances between all pairs
     subplot(8,3,1:3:7); pcolorjk(Mbp_distance); shf; hold on; %plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.');
        axis equal off; set(gca,'ydir','normal','xdir','reverse'); colormap(gca,cm_distance); colorbar('fontsize',12); title('Bipolar pair distance (mm)','fontsize',sizeoffont,'fontweight','normal')
        clim([0 80]);
    % plot confusion matrix of bipolar distances between all pairs
     subplot(8,3,2:3:8); pcolorjk(rad2deg(Mbp_angle)); shf; hold on; %plot([1 1; 1 256]',[1 256;256 256]','k-','linewidth',1); plot([1 1 1 128 256],[1 128 256 256 256],'k.'); %plot([1 128 256 256 256]/2,[1 1 1 128 256]/2,'k.');
        axis equal off; set(gca,'ydir','normal','xdir','reverse'); colormap(gca,hsv(360)); clim([-180 180]); colorbar('fontsize',12); title('Bipolar pair angle (deg)','fontsize',sizeoffont,'fontweight','normal')
    
    % histogram of number of pairs per bin
     subplot(8,2,7); histogram(make1d(Mbp_distance),0:binsz:85,'facecolor',.5*[1 1 1]); set(gca,'fontsize',12); xlabel('Binned bipolar distance (mm)','fontsize',sizeoffont);
     ylabel('Counts/bin','fontweight','normal'); axis tight; grid on; cb=colorbar; set(cb,'visible','off'); xlim(xldist)
    
    
    if save_mat
    % modify path to save_path location
    save(fullfile(data_results, [pt '_fig3_data.mat']), 'pt', 'binz', 'toplot', 'frx', 'binsz', 'Mbp_distance', 'cm_distance');
    end
    
end




