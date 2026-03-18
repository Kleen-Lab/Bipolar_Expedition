function [out] =make1d(in)
% take a whatever-D array and return it as a 1-D vector
% Kleen Lab 2023
out=reshape(in,1,prod(size(in))); %#ok<*PSIZE>

