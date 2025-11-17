function [mz_z,mz,mflat,bpdist ] = bin_zscore_trial(m,nchtocheck,Mbpdist,binsz,frx,none1sqrt2log3)

mflat=[];
bpdist=[];

for c2=1:nchtocheck;
    mflat=[mflat; sq(m(:,c2,:,:))];
    bpdist=[bpdist; Mbpdist(:,c2)]; % corresponding distance index
end

binz=0:binsz:85; % can change binsz PRN at this point to test different bin resolutions on the plots below
%clear mx my mz
nbinz=length(binz)-1;
mz=[];
%bin by distance and take the mean, creating: frequency X distance (binned)
for i=1:nbinz
    mz(:,:,i)=nanmean(mflat(bpdist>=binz(i) & bpdist<binz(i+1),:,:),1);
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

mz_z = [];
mz_zAvg = [];
mz_Avg = squeeze(nanmean(mz,2));
nfrx=length(frx);
% zscore the log transformed power according to frequency
for i=1:nfrx; nns=squeeze(~isnan(nanmean(mz(i,:,:),2)));
    switch none1sqrt2log3
        case 2
            mz_z(i,:,nns)=zscore(sqrt(squeeze(mz(i,:,nns))),[],2);
        case 3
            mz_z(i,:,nns)=zscore( log(squeeze(mz(i,:,nns))),[],2);
    end
end;
% zscore the log transformed power according to frequency
for i=1:nfrx; nns=~isnan(mz_Avg(i,:));
    switch none1sqrt2log3
        case 2
            mz_zAvg(i,nns)=zscore(sqrt(mz_Avg(i,nns)));
        case 3
            mz_zAvg(i,nns)=zscore( log(mz_Avg(i,nns)));
    end
end;

