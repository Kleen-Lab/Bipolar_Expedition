
% BIPOLAR PAIR ANALYSIS: EACH VS. ALL
% Saves data for every patient for specified component: grid, strip, depth
% Data in fix_{component}_tent.mat is relative change from referential
% signal
% Can then plot figure 4 (relative change from referential) using
% fig4_new2025.m

pts = {'EC133', 'EC175', 'EC181', 'EC183', 'EC186', 'EC187', 'EC196', ...
    'EC219', 'EC220', 'EC221', 'EC222', 'EC131', 'EC143', 'EC157', 'EC162', 'EC168'};

for p = 1:length(pts)
    EachVsAll_cleaned2025_fxn(pts{p});
end

function [mDiff, mb_m, mARb_m, Mbp_distance] = EachVsAll_cleaned2025_fxn(pt)

    %if ~exist('pt','var')||isempty(pt); pt='EC175'; end %pt='EC175'; % EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)
    %if ~exist('nchtocheck','var')||isempty(nchtocheck); nchtocheck=128*2; end
    %if ~exist('windowstocheck','var')||isempty(windowstocheck); windowstocheck=250; end %each window is 1 second of data (non-overlapping) %pt='EC175'; % EC175 and EC183 both have intact 16x16 square grids (channel #s 1:256)


    g1s2d3=2; % use either grids (1) or strips (2) or depths (3) but not the others
    doanglerange=0;
    recordings = [];
    none1sqrt2log3=2;
    
    % add paths to helper functions
    
    if g1s2d3 == 1
        component = 'grids';
    elseif g1s2d3 == 2
        component = 'strips';
    else
        component = 'depths';
    
    cm=cool(6); cm(1,:)=[0 0 0]; 
    
    data_root = '/data/';
    datadir = fullfile(data_root, 'baseline-high-density-data');
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
        recordings = [recordings size(nonspks_windows, 2)];
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
    
    an_electrode_info_path = fullfile(data_root, 'bipolar_expedition', 'AN_ElectrodeInfoTDT.xlsx');
    [bpN,bpT]=xlsread(an_electrode_info_path, pts{p});
    [em,eleclabels,anatomy]=getelecs(pts{p},2);
    
    cm=cool(6); cm(1,:)=[0 0 0]; 
    
    datadir = fullfile(datadir, 'bandpassfiltered');
    
    if g1s2d3 == 1
        use_ch = find(strcmpi('grid', anatomy(:,3)) | strcmpi('minigrid', anatomy(:,3)));
        binsz=2; % bin size in mm
        xldist = [0 60];
    elseif g1s2d3 == 2
        use_ch = find(strcmpi('strip', anatomy(:,3)));
        binsz = 4;
        xldist = [0 60];
    else
        use_ch = find(strcmpi('depth', anatomy(:,3)));
        binsz = 4; %binsz=2;
        xldist = [0 60];
    end
    
    badchI=isnan(mean(mean(d,1),3))'; %any channel with nans in any window will be a "bad channel"
    badchI(badchidx{p})=1; %follow up with all marked as bad channels in original preprocessing
    x=[]; xbch=true(size(badchI)); for i=1:size(bpT,1); x=[x bpN(i,1):bpN(i,2)]; end; xbch(x)=false; 
    badchI(xbch)=1; %find and nan empty channels unaccounted for by component rows
    d(:,badchI,:)=nan;
    okc=~badchI; clear x xbch
    
    %windowstocheck=min([windowstocheck size(d,3)]);
    %windowstocheck=1:windowstocheck; %convert to a vector of windows, 1:X
    max_windows = min(250, size(d, 3)); % here is where we specify how many windows
    windowstocheck = 1:max_windows;
    
    % ALL PAIRS (each vs. all others) analysis and example plot
     d=d(:,:,windowstocheck); clear Straces_allch; %free up RAM by getting rid of whatever won't be used (only using first ___ number of windows)
                     % ***opportunity here to select speech or stim windows
     d = d(:, use_ch, :);
     nchtocheck = size(d, 2);
    

     
    
    % euclidean distance and angle (in sagittal plane) for each pair
    
    for c1=1:nchtocheck; 
         for c2=1:nchtocheck; 
             Mbp_distance(c1,c2)=distance3D(em(c1+addons(c1),:),em(c2+addons(c1),:)); 
             %get angle in sagittal plane for each pair (all L-side pts so no need to flip anyone)
              % this is from the inverse tangent of the vertical (superior/inferior) axis difference divided by
              % horizontal (anterior/posterior) axis difference of the two electrodes
              if c1~=c2; Mbp_angle   (c1,c2)=atan2((em(c2,3)-em(c1,3)),(em(c2,2)-em(c1,2))); end
         end; 
     end
        rmv=Mbp_distance==0; %remove spurious zero distance values
        Mbp_distance(rmv)=nan; 
        Mbp_angle   (rmv)=nan; 
        %referential signal denoted as "0" distance for coding/indexing purposes below
        for c1=1:nchtocheck; 
            Mbp_distance(c1,c1)=0; 
            Mbp_angle   (c1,c1)=nan; % angle will still be nans (to avoid confusion with 0 degrees)
        end 
     % remove mirror image in confusion matrix
     for c2=1:nchtocheck; 
         for c1=c2+1:nchtocheck; 
             M(c1,c2,:)=nan; 
             Mbp_distance(c1,c2)=nan; 
             Mbp_angle   (c1,c2)=nan;
         end; 
     end
     Mbp_distance(rmv)=nan; 
     %}
    
    %%
     % Here
     
    %[M,Mrefave,Mbp_distance,frx,~,Mbp_angle]=bpspectra_EachVsAll_2025(d,sfx,frxrange,em,nchtocheck,none1sqrt2log3);
    
    % In main script, before calling the function:
    em_filtered = em(use_ch, :);
    [M,Mrefave,Mbp_distance,frx,~,Mbp_angle]=bpspectra_EachVsAll_2025(d,sfx,frxrange,em_filtered,nchtocheck,none1sqrt2log3);
    
    % [M,Mrefave,frx,~]=dmod_bpspectra_EachVsAll_2023(d,sfx,frxrange,em,nchtocheck);
    
         % M is bipolarchannel1 X bipolarchannel2 X frx X 1secwindow
    nfrx=length(frx);
      
    
    % ANGLE RANGE: if desired, subselect bipolar angle within a given range
    
    if doanglerange
        m=M; Md=Mbp_distance; Ma=Mbp_angle; %make a copy (only do this once) in case you want to repeat later with different angle range (use next chunk of code)
        
        anglemin= 45; anglemax= 90;
        M=m; Mbp_distance=Md; Mbp_angle=Ma; 
        ww=0;
        for c1=1:size(M,1);
        for c2=1:size(M,2);
          if rad2ang(Mbp_angle(c1,c2))<anglemin || rad2ang(Mbp_angle(c1,c2))>anglemax;
            Mbp_distance(c1,c2,:,:)=nan;
            Mbp_angle(c1,c2,:,:)=nan;
            disp([num2str(c1) ' ' num2str(c2)])
            ww=ww+1;
          end
        end
        end
        ww
    end
    
    % Unstack and line up into 2D matrix to create channels^2 X frequencies for easier indexing-->binning
    mb=[]; 
    mARb=[]; 
    binz=[-1 0:binsz:85]; % can change binsz PRN at this point to test different bin resolutions on the plots below
    clear mx my mz  
    nbinz=length(binz)-1;
    
    disp(windowstocheck);
    
    disp(['Size of M: ' mat2str(size(M))]);
    disp(['Size of Mrefave: ' mat2str(size(Mrefave))]);
    disp(['Max window index: ' num2str(max(windowstocheck))]);
    
    for w = windowstocheck
        try
            disp(['Processing window ' num2str(w) '...']);
    
            Mflat = [];
            Mflatbp_distance = [];
            Mrefaveflat = [];
    
            for c2 = 1:nchtocheck
                % Try accessing M and Mrefave safely
                try
                    mc = squeeze(M(:, c2, :, w));
                    mrc = squeeze(Mrefave(:, c2, :, w));
                catch err_inner
                    error(['⚠ERROR accessing M or Mrefave at c2 = ' num2str(c2) ', w = ' num2str(w) ': ' err_inner.message]);
                end
    
                Mflat = [Mflat; mc];
                Mflatbp_distance = [Mflatbp_distance; Mbp_distance(:, c2)];
                Mrefaveflat = [Mrefaveflat; mrc];
            end
    
            % Bin averaging
            for i = 1:nbinz
                mb(:, i, w) = nanmean(Mflat(Mflatbp_distance > binz(i) & Mflatbp_distance <= binz(i+1), :), 1);
                mARb(:, i, w) = nanmean(Mrefaveflat(Mflatbp_distance > binz(i) & Mflatbp_distance <= binz(i+1), :), 1);
            end
    
            disp([num2str(round(w / windowstocheck(end) * 100, 1)) '% done']);
    
        catch err
            disp(['Error at window w = ' num2str(w)]);
            disp(getReport(err, 'extended'));
            return;
        end
    end
    
    disp('Finished all windows');
    
    %{
    parfor w=windowstocheck; disp(num2str(w)); % parfor here to run the windows
        Mflat=[];
          Mflatbp_distance=[];
          %Mflatbp_angle=[];
        Mrefaveflat=[];
        for c2=1:nchtocheck; 
            Mflat=[Mflat; squeeze(M(:,c2,:,w))];               
            Mflatbp_distance=[Mflatbp_distance; Mbp_distance(:,c2)]; % corresponding distance index
            %Mflatbp_angle   =[Mflatbp_angle;    Mbp_angle(:,c2)]; % corresponding angle index
            Mrefaveflat=[Mrefaveflat; squeeze(Mrefave(:,c2,:,w))];              
        end
        
        %bin by distance and take the mean, creating: frequency X binned distance
           % first for referential (distance = 0)
        for i=1:nbinz % binz will be including >lower bound and up to and including (<=) upper bound
            mb  (:,i,w)=nanmean(Mflat      (Mflatbp_distance>binz(i) & Mflatbp_distance<=binz(i+1),:),1); 
            mARb(:,i,w)=nanmean(Mrefaveflat(Mflatbp_distance>binz(i) & Mflatbp_distance<=binz(i+1),:),1);
        end
        
        disp([num2str(round(windowstocheck(w)/windowstocheck(end)*100,1)) '% of windows'])
    end; disp('Done')
    %}
    
    binz(1)=[]; 
    
    %% zscore the log transformed power according to frequency
    
    mb__m=squeeze(mean(sqrt(mb),3)); 
    mARb__m=squeeze(mean(sqrt(mARb),3)); 
     nfrx=length(frx);
     mb__m_z=mb__m; mARb__m_z=mARb__m;
     for i=1:nfrx; nns=~isnan(mb__m(i,:)); 
         mb__m_z(i,nns)=zscore((mb__m(i,nns))); 
         mARb__m_z(i,nns)=zscore((mARb__m(i,nns))); 
     end; 
    

    
    if doanglerange
        subplot(8,6,22); angs=make1d(Mbp_angle); nns=~isnan(angs);
     polarhistogram(angs(nns),-pi:2*pi*(5/360):pi,'facecolor',.5*[1 1 1]); title('Bipolar angle distribution (5^o steps)','fontsize',14,'fontweight','normal'); 
     maxrt=max(get(gca,'rtick')); 
     set(gca,'rtick',[min(get(gca,'rtick')) maxrt/2 maxrt],'ThetaTick',[0 90 180 270],'fontSize',9); 
    end
    
    if doanglerange;  
    % plot illustration of location of bipolar pairs, colored by same distance scale
     subplot(2,3,3); hold on; 
     elecsbrain(pt,0,[1:nchtocheck],[0 0 0],'l',0,5,2); alpha 0.05; litebrain('r',0); zoom(1.5)
     for c1=1:size(Mbp_angle,1)
       for c2=1:size(Mbp_angle,2)
           if ~isnan(Mbp_angle(c1,c2))
                plot3([em(c1,1) em(c2,1)],[em(c1,2) em(c2,2)],[em(c1,3) em(c2,3)],'-','color',cm_distance(1+round(Mbp_distance(c1,c2)),:),'LineWidth',.5)%,'LineWidth',.75)
           end
       end
     end
    end
    
    
    %% bipolar power minus mean referential power for all pairs
    % will also add a line at 10mm for visualization of this common clinical inter-electrode distance
    
    % transform power before calculating percent change % 1: no transform, 2: square root, 3: log
    mb_=mb; mARb_=mARb; % make a copy... then transform the copy
    
    if none1sqrt2log3==1;     txtyp='raw'; % no transform, power in raw form
        mb_(~isnan(mb_))      =    (mb_(~isnan(mb_)));
        mARb_(~isnan(mARb_))  =    (mARb_(~isnan(mARb_)));
    elseif none1sqrt2log3==2; txtyp='square root'; % square root transform
        mb_(~isnan(mb_))      =sqrt(mb_(~isnan(mb_)));
        mARb_(~isnan(mARb_))  =sqrt(mARb_(~isnan(mARb_)));
    elseif none1sqrt2log3==3; txtyp='natural log'; % log transform --> problematic because taking % change of certain small negative values creates extreme values
        mb_(~isnan(mb_))      =log (mb_(~isnan(mb_))+1);
        mARb_(~isnan(mARb_))  =log (mARb_(~isnan(mARb_))+1);
    end
        
    % mean across windows
    mb_m=squeeze(mean(mb_,3)); 
    mARb_m=squeeze(mean(mARb_,3)); 
    mDiff=((mb_m-mARb_m)./mARb_m)*100;
    distance = Mbp_distance;
    
    save(['/Volumes/SPIKE/bipolar_project/bipolar_SeaHorse/devkrish/OMNIRERUN/fix_' component '_tent_' pt '.mat'], 'mDiff', 'mb_m', 'mARb_m');
    save(['/Volumes/SPIKE/bipolar_project/bipolar_SeaHorse/devkrish/OMNIRERUN/fix_' component '_dist_tent_' pt '.mat'], 'distance');

    end

end
