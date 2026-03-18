function pcolorjk(a,b,c)
% adds a row and column of nans so that pcolor can be used w/o losing row
% and column from visualization

if     nargin==1
    a=[a  nan(size(a,1),1)];
    a=[a; nan(1,size(a,2))];
    pcolor(a); shf
elseif nargin==2
    error('Only two inputs to pcolor/pcolorjk, need 1 or 3')
elseif nargin==3
    a=[a a(end)+diff(a(end-1:end))];
    b=[b b(end)+diff(b(end-1:end))];
    c=[c  nan(size(c,1),1)];
    c=[c; nan(1,size(c,2))];
    pcolor(a,b,c); shf
end


