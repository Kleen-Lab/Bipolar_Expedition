function [M,Mrefave,Mbp_distance,frx,Mc,Mbp_angle]=bpspectra_EachVsAll_2025(d,sfx,frxrange,em,nchtocheck,none1sqrt2log3)
% for each window in d (samples x electrodes x windows), take channel c1
% minus channel c2 and computes spectrum on that bipolar signal, placing
% into a confusion matrix. Also does referential spectra (on the diagonal).
% INPUTS:
%   d is a matrix of ICEEG voltage values as samples x electrodes x windows
%   sfx is sampling frequency
%   frxrange is the range of frequencies to look at [min max]
%   em is the XYZ coordinates of the electrodes in d
%   nchtocheck is the number of electrodes (channels) to check if not wanting to do all
% OUTPUTS:
%   M is bipolarchannel1 X bipolarchannel2 X frx X 1secwindow
%   Mbpdist is a confusion matrix of the euclidean distances for every bipolar pair (channel x channel)
%   frx is the frequency index for M (3rd dimension)

%just getting frequency index
 [~,frx]=spectrogramjk_chronuxmtfft(zeros(1,size(d,1)),sfx,frxrange,[.5,1],0); 
 
 nwindtocheck=size(d,3);
 
 % Create matrix M:  channels X channels X frequencies X windows 
 M=nan(nchtocheck,nchtocheck,length(frx),nwindtocheck); 
 Mrefave=M; %copy
 Mbp_distance=nan(nchtocheck,nchtocheck); % will record the euclidean distance for each of these pairs
 Mbp_angle=nan(nchtocheck,nchtocheck); % will record the angle (in sagittal plane relative to horizontal plane as zero degrees) for each of these pairs
 Mc=nan(nchtocheck,length(frx),nwindtocheck); % will record the individual spectrum results
 tic
 for w = 1:nwindtocheck 
     parfor c1=1:nchtocheck % parfor here, seems to cause errors
          % get spectrum of channel 1 (referential)
            fft = spectrogramjk_chronuxmtfft(sq(d(:,c1,w)),sfx,frxrange,[.5,1],0);
            Mc(c1,:,w) = fft;
         for c2=1:nchtocheck 
             % get bipolar trace from channel 1 minus channel 2
             trc=sq(d(:,c1,w))-sq(d(:,c2,w)); 
             if ~any(isnan(trc))
                 % get spectrum of bipolar channel above
                 M(c1,c2,:,w)=spectrogramjk_chronuxmtfft(trc,sfx,frxrange,[.5,1],0); 
                 % take mean of referential spectra from channel 1 and
                 % channel 2 while in log or square root form, then
                 % transform back using exponent or square, respectively
                 if     none1sqrt2log3==1
                   Mrefave(c1,c2,:,w)=...
                          (mean([    (spectrogramjk_chronuxmtfft(sq(d(:,c1,w)),sfx,frxrange,[.5,1],0))...
                                     (spectrogramjk_chronuxmtfft(sq(d(:,c2,w)),sfx,frxrange,[.5,1],0));],   2));
                 elseif     none1sqrt2log3==2
                   Mrefave(c1,c2,:,w)=...
                          (mean([sqrt(spectrogramjk_chronuxmtfft(sq(d(:,c1,w)),sfx,frxrange,[.5,1],0))...
                                 sqrt(spectrogramjk_chronuxmtfft(sq(d(:,c2,w)),sfx,frxrange,[.5,1],0));],   2)).^2;
                 elseif none1sqrt2log3==3
                   Mrefave(c1,c2,:,w)=...
                       exp(mean([ log(spectrogramjk_chronuxmtfft(sq(d(:,c1,w)),sfx,frxrange,[.5,1],0))...
                                  log(spectrogramjk_chronuxmtfft(sq(d(:,c2,w)),sfx,frxrange,[.5,1],0));],   2));
                 end
             end  
         end
     end; disp([num2str(round(w/nwindtocheck*100,1)) '% of windows done, elapsed time: ' num2str(toc)])
 end; clear c1 c2 trc
 % include referential signal power on the diagonal for storage/plotting
 for w = 1:nwindtocheck 
     for c1=1:nchtocheck
             trc=sq(d(:,c1,w));
             if ~any(isnan(trc))
                 M(c1,c1,:,w)=spectrogramjk_chronuxmtfft(trc,sfx,frxrange,[.5,1],0);
             end  
     end;
 end; 

 % euclidean distance and angle (in sagittal plane) for each pair
 for c1=1:nchtocheck; 
     for c2=1:nchtocheck; 
         Mbp_distance(c1,c2)=distance3D(em(c1,:),em(c2,:)); 
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

