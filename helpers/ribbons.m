function [m,envt,envb]=ribbons(x,mtx,c,tr,v,sm,newfig)
% Jon Kleen jon.kleen@ucsf.edu
% will take a matrix of trials (rows) by time (columns) and create
    % envelopes of mean +/- variance (default SEM, can specify SD with v)
%x= The horizontal data points (ie time stamps), can leave empty if desired
%mtx= the matrix of trials (rows) by time (columns)
%c is color (RGB vector)
%tr= transparency, a value ranging from 0 for invisible to 1 for opaque
%v= string specifying SEM or SD or CI for the variance envelope to use
% if size(mtx,1)>size(mtx,2); mtx=mtx'; end
if nargin==1; mtx=x; x=[]; end
if ~exist('x','var') || isempty(x); x=1:size(mtx,2); else x=x(:)'; end
if ~exist('c','var') || isempty(c); c=[0 .75 .75]; end %also nice is [0 .75 1]
if ~exist('tr','var') || isempty(tr); tr=.5; end
if ~exist('v','var') || isempty(v); v='sem'; else v=lower(v); end
if ~exist('sm','var') || isempty(sm); sm=0; elseif ~isnumeric(sm); error('smooth parameter must be an integer'); end
if ~exist('newfig','var') || isempty(newfig); newfig=0; end
m=nanmean(mtx,1); 
envt=nanstd(mtx,[],1); 
if sm; m=smooth(m,sm)'; envt=smooth(envt,sm)'; end
envb=envt;
if strcmpi(v,'sem'); envt=envt/sqrt(size(mtx,1)); envb=envt;
elseif strcmpi(v,'ci'); envt=prctile(mtx,97.5)-m; envb=m-prctile(mtx,2.5); 
end
edges=c-.2; edges(edges<.1)=.1; mc=c-.2; mc(mc<0)=0;
if newfig; figure; end
fill([x fliplr(x)],[m+envt fliplr(m-envb)],c,'EdgeColor',edges,'FaceAlpha',tr,'EdgeAlpha',tr);
if ~isempty(m); hold on; plot(x,m,'color',mc,'linewidth',.25); hold off; end
