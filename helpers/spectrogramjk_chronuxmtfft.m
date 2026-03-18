function [S,frx]=spectrogramjk_chronuxmtfft(data,sfx,frxrange,movingwindow,fig,logspacefrx)

% Need chronux toolbox installed: http://chronux.org/
% 
% Example:
% figure('color','w')
% sp(5,1,1:4); spectrogramjk_chronuxmtfft(data,sfx,[0 255],0)
% sp(5,1,5); plot(data); axis tight


if ~exist('frxrange','var'); frxrange=[0 sfx/2]; end
if ~exist('fig','var'); fig=1; end 
if ~exist('logspacefrx','var'); logspacefrx=0; end 

params.fpass=frxrange;
params.tapers=[3 5];
params.Fs=sfx;

if ~exist('movingwindow','var')
winsize=1/4; %in seconds
winstep=1/16; %in seconds
movingwindow=[winsize winstep];
end

if logspacefrx % not working well yet, some frequencies empty
    %[cfrx,~]=filterbankjk(frxrange,1/7); 
    cfrx=[2,3,4,5,6,7,8,10,12,15,18,22,27,33,41,50,60,74,90,110,134,163,198];
    for i=1:length(cfrx-1); params.fpass=cfrx([i i+1]);
        [s,~,~]=mtspecgramc(data,movingwindow,params);
        S(i)=mean(s);
    end
    frx=cfrx;
else
    [S,t,frx]=mtspecgramc(data,movingwindow,params);
    S=S';
end



if fig; figure('color','w'); 
pcolor(t,frx,log10(S)); shf
title('log10(power)')
xlabel('Time (s)'); ylabel('Frequency (Hz)');
end
% subplot(121)
% pcolor(t,f,log(S')); shf
% subplot(122)
% plot_matrix(S,t,f);
% xlabel([]); % plot spectrogram
% caxis([8 28]); colorbar;