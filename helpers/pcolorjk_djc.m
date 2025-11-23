function pcolorjk_djc(x,y,z)
% adds a row and column of nans so that pcolor can be used w/o losing row
% and column from visualization


x=[x  nan(size(x,1),1)];

y=[y  nan(size(y,1),1)];

z=[z  nan(size(z,1),1)];
z=[z; nan(1,size(z,2))];


pcolor(x,y,z); shf


