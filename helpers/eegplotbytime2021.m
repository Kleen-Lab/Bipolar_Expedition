function []=eegplotbytime2021(d,sfx,ampfactor,strt,fig,colr,lw)
% displays EEG data channels by time
% jon.kleen@ucsf.edu

% if size(d,1)>size(d,2); d=d'; end
[ch,t]=size(d);
% checking for inputs, making defaults if user didn't include
if ~exist('ampfactor','var'); ampfactor=3; end % gain (will multiply by this number)
    ampfactor(ampfactor<1)=1;
if ~exist('strt','var')||isempty(strt); strt=0; end %jump ahead to start time "strt" in seconds
if exist('fig','var') && fig; figure; end
if ~exist('colr','var'); colr='b'; end
if ~exist('lw','var'); lw=.25; end

% center each lead at zero
for i=1:size(d,1);
    pc=prctile(d(i,:),50); % perc is percentile to amplify by
    d(i,:)=d(i,:)-pc;
end

m1d = reshape(d,1,size(d,2)*size(d,1));
maxd=max(m1d);
am=((ampfactor/10000)*maxd)/maxd;
shift = repmat(-1*(1:ch)',1,t); 

x=0:1/sfx:t/sfx-1/sfx; x=x+strt; %time vector, shift to start time

myplot=plot(x,d*am+shift,'color',colr,'linewidth',lw);
ylim([-1*size(d,1)-1 0])
axis tight; xlabel('time (sec)')
hold off
