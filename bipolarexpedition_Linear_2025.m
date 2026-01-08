% BIPOLAR PAIR ANALYSIS: LINEAR
% Display linear bipolar analysis plots and linear analysis


% Section 1: linear_analysis
% Run linear analysis on grids, strips, or depths individually, and save
% the spectral data per distance into separate arrays

%function linear_analysis(component_num,transform)
%component_num --> grids(1), strips(2), depths(3)
%transform --> use 2 for sqrt


for q = 1:3 %run 1:3 for all components
    datadir = getenv("BIPOLAR_DATA");
    
    save_data = true; % save plots averaged across all patients for each component spectra
    save_data_pt = true; % save plots for example patient spectra
    save_data_path = fullfile(datadir,'results'); %update to path of output
    component_num = q;
    transform = 2;
    g1s2d3=component_num;
    %extra_plots = false;
    none1sqrt2log3 = transform;


    % how many bipolar steps discretized distances to loop through
    if g1s2d3==1;     maxbpd=5;
    elseif g1s2d3==2; maxbpd=2;
    elseif g1s2d3==3; maxbpd=4; %change to 4
    end

    switch g1s2d3
        case 1; bpd_mm=4;
        case 2; bpd_mm=10;
        case 3; bpd_mm=5;
    end
    bpd_mm=bpd_mm*(0:10); %bipolar distances to be evaluated, in mm

    caxisrange=[0 20];
    cm=cool(17);
    cm=[0 0 0;1 1 1;1 1 1;cm]; %first entry black for referential, rest allows color-coding of physical distance

    tag_spikes_path = fullfile(datadir, 'taggedspikes_April2022.mat');
    disp('loading tagged');
    load(tag_spikes_path);
    sfx=512;
    frxrange=[2 200]; %frequency range to examine
    ft=[2 5 10 20 50 100 200]; %frequency tick marks for plots
    ftl=cellstr(num2str(ft')); %frequency labels for plots

    %which patients are ok to do
    okpt=false(1,length(pts));
    okpt([3 4 8 9 10 11 12:16 19:23])=1;
    newpts = pts(okpt);

    %% LINEAR PAIRS analysis and plots
    figure(1); set(gcf,'color','w','position',[372 1 1297 1337]);
    u=dir(fullfile(datadir, 'baseline-high-density-data')); uptbl={};
    for i=1:length(u)
        uname=u(i).name;
        uptbl{i,1}=uname(1:end-28);
    end

    uptbl(1:2)=[];
    clear i u uname
    hasmat=false(maxbpd+1,length(pts));

    trm_data = struct();
    trm_data.grids = cell(6, length(pts));   % 6 distances: 0,4,8,12,16,20 mm
    trm_data.strips = cell(3, length(pts));  % 3 distances: 0,10,20 mm
    trm_data.depths = cell(5, length(pts));  % 5 distances: 0,5,10,15,20 mm
    trm_data.distances_mm.grids = [0, 4, 8, 12, 16, 20];
    trm_data.distances_mm.strips = [0, 10, 20];
    trm_data.distances_mm.depths = [0, 5, 10, 15, 20];
    trm_data.frequency = []; % Will be filled with frx
    trm_data.patients = pts;

    for p=find(okpt) %[4 12:23]
        pblocks=strfind(uptbl,pts{p});
        for i=1:length(pblocks)
            isbl(i,1)=~isempty(pblocks{i});
        end
        ptbl=find(isbl);
        if ~isempty(ptbl); disp(['Loading ' pts{p} ' blocks...']); end

        d=[]; %d will become a matrix of samples by channels by trials consisting of referential intracranial EEG data
        nwind=0;
        for b=1:length(ptbl); disp(uptbl{ptbl(b)})
            datapath = fullfile(datadir, 'baseline-high-density-data', [uptbl{ptbl(b)} '_baselineWindows_fromraw.mat']);
            load(datapath);
            % get rid of baseline windows containing spikes or artifact
            spksarti=hasspk | hasarti;
            nonspks_windows(:,spksarti)=[];
            hasstim(spksarti)=[];
            hasspeech(spksarti)=[];
            clear hasspkvec hasspk hasartivec hasarti spksarti % now clear spike- and artifact-related variables from workspace

            % convert to 3D matrix, combine all windows from consecutive blocks for each patient
            for i=1:size(nonspks_windows,2)
                d(:,:,i+nwind)=nonspks_windows{2,i}'; % d is a 3D matrix samples by channels by trials
            end
            nwind=size(d,3);

            clear nonspks_windows info
        end; clear b

        nch=size(d,2);


        %% get electrode info
        % bpN is the first (column 1) and last (column 2) electrodes for
        %   each component (grid strip or depth). Column 3 is for grids only,
        %   gives the number of electrodes in each row.
        % bpT is the name of the component (column 1) and whether it is a grid strip or depth
        %[bpN,bpT]=xlsread(['/Users/davidcaldwell/code/high_density_ecog/AN_ElectrodeInfoTDT.xlsx'],pts{p});
        an_electrode_info_path = fullfile(datadir, 'AN_ElectrodeInfoTDT.xlsx');
        [bpN,bpT]=xlsread(an_electrode_info_path,pts{p});

        % pull electrode XYZ coordinates (em) from TDT_elecs_all.mat file
        [em,~,~]=getelecs(pts{p},2);
        if isempty(em); error(['Missing electrode coordinates for ' pts{p} ', check paths']); end

        D=d;


        for bpd=0:maxbpd %bipolar distance (# of electrodes to subsample)
            %% bipolar conversion
            % d is a matrix of samples by channels by trials
            %   initially it consists of referential intracranial EEG data
            %   then below it gets converted to bipolar re-referenced data,
            %   with nan values replacing the entries where bipolar pairs went past the end
            %   for example, a usual bipolar subtracts each electrode from the
            %   next and so the last electrode in the component/row will be nan
            %   (eg., 8 ref elecs --> skip 0 bp --> 7 bp elecs + nan entry)
            %   (eg., 8 ref elecs --> skip 1 bp --> 6 bp elecs + 2 nan entries), etc

            d=D; %REFRESH with original referential data

            bp_distance=nan(nch,1); %will fill in euclidean distance in 3D space for each bipolar pair created
            if bpd>0
                % linear grid rows, strips, depths
                for jj = 1:size(d,3) %windows
                    for r=1:size(bpT,1) %each row of the sheet is a component (grid, strip, or depth)
                        if any(strcmpi(bpT(r,2),{'grid','minigrid'})) %grids (2-D)
                            N=bpN(r,3)-bpd; % Number of new consecutive bipolar contacts of distance "bpd"
                            if N>0
                                for i=bpN(r,1):N+bpd:bpN(r,2); %every grid row
                                    c1=[i:i+N-1];
                                    c2=[i:i+N-1]+bpd;
                                    d(:,c1,jj)=d(:,c1,jj)-d(:,c2,jj);
                                    d(:,i+N:i+N+bpd-1,jj)=NaN; %last channel in the line will be NaNs
                                    for i=1:length(c1);
                                        % get actual distance in 3D space
                                        bp_distance(c1(i))=distance3D(em(c1(i),:),em(c2(i),:));
                                        %get angle in sagittal plane (all L-side pts so no need to flip anyone)
                                        bp_angle   (c1(i))=atan((em(c2(i),3)-em(c1(i),3)) / (em(c2(i),2)-em(c1(i),2)));
                                    end
                                end
                            else %if bipolar spacing is longer than #electrodes in component, nan the component electrodes
                                d(:,bpN(r,1):bpN(r,2),:)=nan;
                                bp_distance(:,bpN(r,1):bpN(r,2))=nan;
                            end
                        else; i=bpN(r,1); %strips, depths (1-D)
                            N=diff(bpN(r,1:2))+1-bpd; % number of new consecutive bipolar contacts of distance "bpd"
                            if N>0
                                c1=[i:i+N-1];
                                c2=[i:i+N-1]+bpd;
                                d(:,c1,jj)=d(:,c1,jj)-d(:,c2,jj);
                                d(:,i+N:i+N+bpd-1,jj)=NaN; %last channel in the line will be NaNs
                                for i=1:length(c1); bp_distance(c1(i))=distance3D(em(c1(i),:),em(c2(i),:)); end
                            else %if bipolar spacing is longer than #electrodes in component, nan the component electrodes
                                d(:,bpN(r,1):bpN(r,2),:)=nan;
                                bp_distance(bpN(r,1):bpN(r,2))  =nan;
                                bp_angle   (bpN(r,1):bpN(r,2))  =nan;
                            end
                        end; clear N i
                    end
                end; clear r jj
            end

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


            %% Calculate spectra and put into matrices (bipolarDistance X patient X frequency) aggregated for each patient
            if any(okc)
                [s,frx]=bpspectra_Linear_2025(d,sfx,frxrange,okc);

                % transform power before analyzing
                s_Tx=s; % make a copy... then transform the copy
                if none1sqrt2log3==1;     txtyp='raw'; % no transform, power in raw form

                elseif none1sqrt2log3==2; txtyp='sqrt'; % square root transform
                    s_Tx=sqrt(s_Tx);
                elseif none1sqrt2log3==3; txtyp='natural log'; % log transform --> problematic because taking % change of certain small negative values creates extreme values
                    s_Tx=log(s_Tx);
                end

                trm=squeeze(mean(s_Tx,3))'; %"trial mean": trm is mean of of power across windows/trials
                trm(~okc,:)=nan;

                if g1s2d3 == 1  % grids
                    trm_data.grids{bpd+1, p} = trm(okc,:);
                elseif g1s2d3 == 2  % strips
                    trm_data.strips{bpd+1, p} = trm(okc,:);
                elseif g1s2d3 == 3  % depths
                    trm_data.depths{bpd+1, p} = trm(okc,:);
                end

                if isempty(trm_data.frequency)
                    trm_data.frequency = frx;
                end

                TRM(bpd+1,p,:)=mean(trm(okc,:),1); % Mean
                TRSD(bpd+1,p,:)=std(trm(okc,:),[],1); %Standard deviation
                TRSE(bpd+1,p,:)=TRSD(bpd+1,p,:)/sqrt(sum(okc)); %SEM
                TRbp_distance{bpd+1,p}=bp_distance; %actual euclidean distances for each bp pair

                SS{bpd+1,p}=s_Tx;
                hasmat(bpd+1,p)=1;
                Sokc{bpd+1,p}=okc;

                figure(2 + double(p==4)); set(gcf,'position',[372 1 1297 865],'color','w');
                if p~=4; sp(5,5,p); else; set(gcf,'position',[495 308 442 453]); end
                hold on; %each patient in their own plot



                % THIS IS WHERE YOU PUT RIBBONS DEVON DAMN IT
                ribbons(frx,trm(okc,:),cm(max([1 bpd_mm(bpd+1)]),:),.5,'sem',0,0); set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl)
                grid on; title(pts{p}); drawnow;
                if p==4;
                    xlabel('Frequency (Hz)');
                    ylabel([txtyp '(power)']);
                    legend({'referential','', num2str(bpd_mm(2)),'', num2str(bpd_mm(3)),'',num2str(bpd_mm(4)),'',num2str(bpd_mm(5)),'',num2str(bpd_mm(6))},'location','sw');
                    axis tight; set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl)
                    colormap(cm);
                    caxis(caxisrange);
                    cb=colorbar;
                    cb.Ticks=[0.5 bpd_mm(2:end)-.5];
                    cb.TickLabels=[{'Referential'} ; cellstr(num2str(bpd_mm(2:end)'))];
                end
            end

            %% ECoG trace plots for increasing bipolar spacing (example patient)
            if p==4
                figure(5); set(gcf,'color','w','position',[1 187 1586 860])
                %sp(6,2,(bpd+1)*2-1);
                sp(4,6,bpd+1);
                chtoplot=49:64; %example channels [EC143-49:64, EC175-]
                windowtoplot=25; %example window
                ts=1/sfx:1/sfx:1; %timestamps for 1-sec window
                eegplotbytime2021(d(:,chtoplot,windowtoplot)',sfx,250,[],0,[.3 .3 .3],1);
                %             for c=1:length(chtoplot); plot(ts,-c*1000+d(:,chtoplot(c),windowtoplot),'color',[0 .6 .6],'linewidth',1); end
                if ~exist('yl','var'); yl=ylim; end; ylim(yl);
                ylim(yl+(bpd/2)) %ylim(yl+(bpd/2+.5))
                axis off
                sp(2,3,4); hold on;
                for c=chtoplot
                    [specttoplot(c-chtoplot(1)+1,:),frx]=spectrogramjk_chronuxmtfft(squeeze(d(:,c,windowtoplot)),sfx,frxrange,[.5,1],0);
                end

                % transform power before plotting
                if none1sqrt2log3==1
                    specttoplot=    (specttoplot);  % no transform, power in raw form
                elseif none1sqrt2log3==2 % square root transform
                    specttoplot=sqrt(specttoplot);
                elseif none1sqrt2log3==3 % log transform --> problematic because taking % change of certain small negative values creates extreme values
                    specttoplot= log(specttoplot);
                end

                if (g1s2d3 == 1) & (save_data_pt)
                    if bpd == 0
                        dk_spec = cell(1, maxbpd + 1);  % Initialize on first distance
                    end
                    dk_spec{bpd+1} = specttoplot;
                end

                ribbons(frx,specttoplot,cm(max([1 bpd_mm(bpd+1)]),:),.3,'sem',0,0);

                %ribbons(frx,trm(chtoplot,:),cm(max([1 bpd_mm(bpd+1)]),:),.3,'sem',0,0);
                grid on;
                set(gca,'xlim',frxrange,'xscale','log','xtick',ft,'XTickLabel',ftl)
                clear c specttoplot
                ylabel([txtyp '(power)']); drawnow;

                if bpd==0
                    subplot(2,100,146:175)
                    %elecsbrain(pts{p},0,[],[0 0 0],'l',0,10,2); alpha 1;
                    %elecsbrain(pts{p},0,[chtoplot],[1 0 0],'l',0,7,2); litebrain('l',.2)
                    zoom(1)
                end

                subplot(2,100,176:200)
                colormap(gca,cm);
                caxis(caxisrange);
                cb=colorbar;
                cb.Ticks=[0.5 bpd_mm(2:end)-.5];
                cb.TickLabels=[{'Referential'} ; cellstr(num2str(bpd_mm(2:end)'))];

                axis off

                if save_data_pt
                    save('spec_sqrt.mat', 'dk_spec', 'frx', 'ft', 'ftl', 'frxrange', 'txtyp');
                end

            end %set(gcf,'position',[1198 785 498 481]) %to resize spectra for figure


        end % patient loop
        clear d isbl ptbl pblocks s trm badchI okc nch nwind

   

    end % bipolar spacing loop

    % additional plots for Linear bipolar analysis: across all patients
    % fill nans for empty (non-analyzed) patients
    TRM(:,~okpt,:)=nan;
    TRSD(:,~okpt,:)=nan;
    TRSE(:,~okpt,:)=nan;
    % TRbpdist(:,~okpt)=nan;  %********

    % plots aggregated across patients
    figure(4); set(gcf,'color','w','position',[88 122 494 624]); %,'position',[1055 216 1217 826]
    rebase=1;
    rebase_fl=[2 10]; %frequency limits for rebasing to referential signal
    p_rebased=nan(length(pts),1);
    hold on
    ps=nan(maxbpd,length(pts),length(frx));
    C=nan(maxbpd,length(frx)); Cp=C;

    for p=1:length(pts);
        for bpd=[0 1:maxbpd]
            if hasmat(bpd+1,p);
                % for this current [bpd+1] distance only!
                p_SS_allwindows=SS{bpd+1,p}; % frequencies X electrodes X windows
                p_SSfrxelecs=mean(p_SS_allwindows(:,Sokc{bpd+1,p},:),3); % frequencies by OKelectrodes
                p_SSfrx=mean(p_SSfrxelecs,2); % frequencies
                ps(bpd+1,p,:)=mean(p_SSfrx,2); % distance X patient X frequencies
            end
        end
        if any(hasmat(bpd+1,p))
            % Rebasing to MEDIAN of referential spectrum
            if rebase
                p_rebased(p)=median(ps(1,p, frx>=rebase_fl(1) & frx<=rebase_fl(2)));
                for bpd=[0 1:maxbpd]
                    ps(bpd+1,p,:)=ps(bpd+1,p,:)-p_rebased(p);
                end
            end
        end
    end

    % Saving power spectra averaged across all patients

    %{
    dk_lin = cell(maxbpd+1, 1);
    
    for i = 1:maxbpd+1
        valid_patients = hasmat(i,:);
        if any(valid_patients)
            dk_lin{i} = squeeze(nanmean(ps(i,valid_patients,:), 2)); %100 F x 1`
        else
            dk_lin{i} = nan(size(frx)); % fill with NaNs if no data
        end
    
    end
    %}

    dk_lin = cell(1, maxbpd + 1);

    for i = 1:maxbpd + 1
        dk_lin{i} = squeeze(ps(i, :, :));  % 11 patients × 100 frequency bins
    end

    % ----

    if save_data

        if g1s2d3 == 1
            save([save_data_path 'grid_trm_data.mat'], 'trm_data')
        elseif g1s2d3 == 2
            save([save_data_path 'strip_trm_data.mat'], 'trm_data')
        else
            save([save_data_path 'depth_trm_data.mat'], 'trm_data')
        end

        disp(['Saved trm_data to: ' save_data_path]);


    end

end


for bpd=[0 1:maxbpd]
    if any(hasmat(bpd+1,:))
        %ribbons(frx,sq(ps(bpd+1,find(hasmat(bpd+1,:)),:)),cm(max([1 bpd_mm(bpd+1)]),:),.5,'sem',0,0);
        %plot(frx,sq(nanmean(ps(bpd+1,hasmat(bpd+1,:),:),2)), 'color',cm(max([1 bpd_mm(bpd+1)]),:),'linewidth',4);
        plot(frx,squeeze((ps(bpd+1,hasmat(bpd+1,:),:))), 'color',cm(max([1 bpd_mm(bpd+1)]),:),'linewidth',1);
        hold on
        if bpd>0
            [clusters, p_values, ~, ~ ] = permutest(squeeze(ps(bpd+1,hasmat(bpd+1,:),:))',squeeze(ps(1,hasmat(bpd+1,:),:))',0); clusters(p_values>=0.05)=[]; p_values(p_values>=0.05)=[];
            for i=1:length(clusters); Cp(bpd,clusters{i})=p_values(i); C(bpd,clusters{i})=1; end;
        end
    end
end
grid on; ylabel([txtyp '(power)']); xlabel('Frequency (Hz)'); ttl=['Mean across patients (n=' num2str(length(find(sum(hasmat,1)))) ')']; title({ttl,''});
if rebase; title({ttl,['Rebased to ' num2str(rebase_fl(1)) '-' num2str(rebase_fl(2)) ' Hz referential signal'],''}); end
set(gca,'xscale','log','xlim',frxrange,'xtick',ft,'XTickLabel',ftl,'ylim',[min(ylim)-.5 max(ylim)])
for bpd=1:maxbpd; plot(frx,C(bpd,:)*.05*bpd+min(ylim)+.01,'-', 'color',cm(max([1 bpd_mm(bpd+1)]),:),'linewidth',2); end
colormap(gca,cm);
caxis(caxisrange);
cb=colorbar;
cb.Ticks=[0.5 bpd_mm(2:end)-.5];
cb.TickLabels=[{'Referential'} ; cellstr(num2str(bpd_mm(2:end)'))];

