function [L]=line_length_alg(d,sfx,llw)
%originally by jon.kleen@ucsf.edu 2016-2021
% modified by DJC 2021
%
% Transforms data into linelength then detects events (spikes) surpassing
% the designated percentile threshold. Note that this function assumes any
% detections in any channel occurring simultaneously are involved in the
% same spike event.
% Based on Estellar et al 2001, DOI 10.1109/IEMBS.2001.1020545
%INPUTS
% d: vector or matrix of ICEEG data, 2d -  (time x windows/channels) or 3d - (time x windows x channels), time
% should be first dimension
% sfx: sampling frequency
% llw: linelength window (in seconds) over which to calculate transform
% OUTPUTS
% transformed data
% L: transformed data
% ets: matrix of events (rows) and their on/off times (2 columns) in samples
% ech: logical index of which channels are involved in each detection,
%        thus having the same number of rows (spikes) as ets

%Example: [ets,ech]=LLspikedetector(d,512,.04,99.99)

if ~exist('llw','var'); llw=.04; end %default linelength window for transform is 40ms
if length(size(d))>3; error('Accepts only vector, 2-D, 3-D matrix for data'); end



%%  1. LINE-LENGTH TRANSFORM
% Will be same size as d with tail end (<llw) padded with NaNs.
numsamples=round(llw*sfx); % number of samples in the transform window
if any(size(d)==1)   %if d is a vector
    L=nan(1,length(d)); % will fill this with transformed data in loop below
    for i=1:length(d)-numsamples
        L(i)=sum(abs(diff(d(i:i+numsamples-1))));
    end
elseif numel(size(d)) == 2           %if d is a 2-D matrix
    L=nan(size(d));
    for i=1:size(d,1)-numsamples
        L(i,:)=sum(abs(diff(d(i:i+numsamples-1,:),1,1)),1);
    end
elseif  numel(size(d)) == 3% if d is a 3-D matrix
    L=nan(size(d));
    for i = 1:size(d,1)-numsamples
        L(i,:,:)=sum(abs(diff(d(i:i+numsamples-1,:,:),1,1)),1);
    end


    % Tada
end