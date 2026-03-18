function [tD] = abs_deriv(d)
% take the absolute derivative of a signal along the first dimension of
% 1D, 2D, or 3D matrix
%
%
% David Caldwell 1.4.2022

%%  1. LINE-LENGTH TRANSFORM
% Will be same size as d with tail end (<llw) padded with NaNs.
if any(size(d)==1)   %if d is a vector
    tD=abs(diff(d));
elseif numel(size(d)) >= 2           %if d is a 2-D or 3-Dmatrix
    tD=abs(diff(d,1,1));
end

