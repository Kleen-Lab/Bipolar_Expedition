function [f] = jkfilt(d,sfx,lo,hi,type,fo) % Jon Kleen 2016

% Band-pass filter (for EEG/ECoG data), e.g. f=jkfilt(d,512,4,12,'butter')
% d is data as a vector or 2-d matrix (it assumes longer dimension is time)
% sfx is sampling frequency
% lo and hi are the respective low and high frequency cutoffs for the band
% type is a string that specifies which of 2 filters to use: 
%     butterfield ('butter' or 'butterfield') or cheby type II ('cheby' or 'cheby2')
if length(size(d))>2; disp('too many dimensions. vector or 2-d matrix only'); return; end
flip=0; if size(d,1)<size(d,2); d=d'; flip=1; end

rp=20; %ripple parameter
freqs=[lo hi]/(sfx/2); %frequency parameters

if nargin<6||isempty(fo); fo=4; end %filter order parameter

if nargin<5; type='butter'; 
elseif ~strcmp(type,'cheby') && ~strcmp(type,'cheby2') && ~strcmp(type,'butter'); 
    type='butter'; disp('FYI: defaulting to butterfield')
end; 
if size(d,1)>1 && size(d,2)>1; %if a matrix
    f=zeros(size(d));
        if strcmp(type,'butter')||strcmp(type,'butterfield'); [b,a]=butter(2,freqs,'bandpass'); 
        elseif strcmp(type,'cheby2')||strcmp(type,'cheby'); [b,a]=cheby2(fo,rp,freqs);
        end
    for i=1:size(d,2); %I timed this loop, actually 4x faster than giving it a matrix
        f(:,i)=filtfilt(b,a,d(:,i));
    end
else %if a vector
    if strcmp(type,'butter')||strcmp(type,'butterfield'); [b,a]=butter(2,freqs,'bandpass'); 
    elseif strcmp(type,'cheby2')||strcmp(type,'cheby'); [b,a]=cheby2(fo,rp,freqs);
    end
    f = filtfilt(b,a,d);
end
if flip; f=f'; end %flip processed data to original orientation of d
