function [s,frx]=bpspectra_Linear_2025(d,sfx,frxrange,okc)
% INPUTS:
%   d is a matrix of ICEEG voltage values as samples x electrodes x windows
%   sfx is sampling frequency
%   frxrange is the range of frequencies to look at [min max]
%   okc is an index of the ok channels (electrodes). 1 is ok, 0 is bad channel
% OUTPUTS
%    ln_s is all spectra data (frx X channel X window) that was log-transformed. Can do mean to create trm.
%    frx is index of frequencies

[~,nch,nwind]=size(d);
%just getting frequency index
[~,frx]=spectrogramjk_chronuxmtfft(squeeze(d(:,1,1)),sfx,frxrange,[.5,1],0); 

% FFT
D=d(:,okc,:);
nfrx=length(frx); %number of anticipated frequencies
spect=nan(nfrx,size(D,2),size(D,3));
L=size(D,2);
tic
parfor w=1:nwind
    W=squeeze(D(:,:,w))'; 
    for c=1:L %c=find(~badchAll(b,:))
        spect(:,c,w)=spectrogramjk_chronuxmtfft(W(c,:),sfx,frxrange,[.5,1],0);
    end
end; clear D W c w L
toc
s=nan(nfrx,nch,nwind);
s(:,okc,:)=spect; clear spect 
