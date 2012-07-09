function out = my_smooth(in,n1,n2)
% 2-dimensional smoothing by convolution with triangular filter of length n
%
% use:
%   out = my_smooth(in,n1,n2)
% 
% input:
%   in - matrix.
%   n1 - smoothing length in the first direction (rows)
%   n2 - smoothing length in the second direction (columns)

t1 = [1:n1 n1-1:-1:1]/n1^2;
t2 = [1:n2 n2-1:-1:1]/n2^2;


out = conv2(t1,t2,in,'same');