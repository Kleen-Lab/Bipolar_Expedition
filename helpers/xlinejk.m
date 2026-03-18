function xlinejk(y,spec,lw)
% plot horizontal line(s) emanating from specific points ("y") on the y axis 
% that are overlaid on the current plot with various properties. Examples:
% spec is k- for black solid line. G actually gives grey! (g already taken for green)
% lw is a number representing linewidth)
% Kleen Lab UCSF
if ischar(y); spec=y; clear y; y=[]; end %in case simply a colored line in the middle is wanted
if ~exist('y','var') || isempty(y); ys=ylim; y=mean(ys); end
if ~exist('spec','var') || isempty(spec); spec='k:';end
if         ~ischar(spec) && all(size(spec)==[1 3]); c=spec; spec=':'; 
    elseif length(spec)==2; c=spec(1); spec=spec(2); 
    elseif length(spec)==1; if ~ischar(spec); c='k'; else c=spec; end
    else c='k';
end
    if strcmp(c,'G'); c=[.75 .75 .75]; end
if nargin<3;lw=1;end
x=xlim; hold on; 
for i=1:length(y); 
plot(x,[y(i) y(i)],spec,'linewidth',lw,'color',c)
end
