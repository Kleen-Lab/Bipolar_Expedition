function elecsbrain(pt,numlab,elecs,colr,plotbrain,newfig,msize,clin1tdt2MNI3,shiftoff)
% Plots specified (or all, if unspecified) TDT electrodes on brain or in space
% Inputs: 
%   pt: patient as a string, such as 'EC133'
%   numlab: labels to plot next to electrodes, 0 for no labels, 1 for number labels, 'a' for actual montage labels
%   elecs: enter vector of electrode numbers to plot. To plot all electrodes, leave empty []
%         Example, [1 3 10] plots electrodes # 1, 3, and 10 from the montage designated in "clin1tdt2MNI3" below 
%   colr: color for electrode dots, as RBG values, such as [0 0 1] for blue
%   plotbrain: hemisphere reconstruction plotting. Empty [] for none, 'l' for left hemisphere, 'r' for right, and 'b' for both/bilateral
%   newfig: 1 to make a new figure, 0 to plot in current figure/axis
%   msize: marker size, default 5
%   clin1tdt2MNI3: which electrode file to use. 1 for clinical montage, 2 for TDT, 3 for MNI-based warped locations (freesurfer)
% Example usage: elecsbrain('EC133','a',[],[0 0 0],'b',1,5,2);
%                elecsbrain('EC133','a',[65:128],[0 0 0],'b',1,5,2);
%                elecsbrain('EC133','a',getregionelecs('EC133','mtg'),[0 0 0],'b',1,5,2);
% See also: getbrain.m, getelecs.m, getregionelecs.m and getregionelecs_verified.m


if iscell(elecs)||ischar(elecs); elecs=getregionelecs(pt,elecs); end
if ~exist('clin1tdt2MNI3','var')||strcmp(clin1tdt2MNI3,'c'); clin1tdt2MNI3=1; end
if ~exist('shiftoff','var'); shiftoff=0; end
% [em,el,an]=getTDTelecs(pt,clin1tdt2==1);%-clin1tdt2+1); 
% if isempty(em); [em,el,an]=getOPSCEAelecs(pt,1); end
[em,~,an]=getelecs(pt,clin1tdt2MNI3); if isempty(em); disp('Cant find electrode file, check patient, path and filename, and whether *****_elecs_all.mat file exists in that path'); end

if nargin<2; numlab=1; end %to display electrode number next to each
if nargin<3 || (exist('elecs','var') && isempty(elecs)); elecs=1:size(em,1); end
elecs=unique(elecs);
if nargin<4 || (exist('colr','var') && isempty(colr)); colr='b'; end %specify electrodes' color, blue default
if nargin<5 || (exist('plotbrain','var') && isempty(plotbrain)); plotbrain='b'; end %to plot brain surface too
    if plotbrain~=0 && ~any(strcmp(plotbrain,{'r','l','b'})); error('Please specify r, l, or b for plotbrain argument'); end
if nargin<7 || (exist('msize','var') && isempty(msize)); msize=5; end %to plot brain surface too

if newfig; figure('color','w','Position',[128 163 792 535]); end %exist('plotbrain','var') && plotbrain>0 && 
% if plotbrain>0; getbrain(pt,1,0,plotbrain); %any(strcmp(pt,{'EC132','EC161','EC162','EC171','EC172','EC173'}))
    if ischar(plotbrain)&&~isempty(plotbrain); getbrain(pt,1,0,plotbrain); end
% end

if shiftoff
    em(em(:,1)<0,1)=em(em(:,1)<0,1)-2; %shift L-sided electrodes 2mm to L direction to help visualize for lateral views
    em(em(:,1)>0,1)=em(em(:,1)>0,1)+2; %shift R-sided electrodes 2mm to R direction to help visualize for lateral views
    em(:,3)=em(:,3)-2; %shift off brain surface for inferior views
    % 
    % em(:,1)=em(:,1)*.99; %shift off brain surface for lateral views
    % em(:,2)=em(:,2)+2; %shift off brain surface for inferior views
end

if numlab~=0; if ischar(numlab); labl=an(:,1); else labl=1:size(em,1); end; 
        el_add(em(elecs,:),'color',colr,'msize',msize,'numbers',labl(elecs));
else    el_add(em(elecs,:),'color',colr,'msize',msize,'edgecol',[.2 .2 .2])
end

angl=nanmean(em(:,1));
if angl<-10; view(270, 0); %assumes L-side electrodes
elseif angl>10; view(90, 0); %assumes R-side electrodes
else abs(angl)<10; view(0, 0); %assumes Bilateral electrodes
end

if plotbrain && newfig; alpha(.65); end

set(gca,'Clipping','off')
cameratoolbar('setmode',''); 