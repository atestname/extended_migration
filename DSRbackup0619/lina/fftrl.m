function [b,f] = fftrl(a,t,mode)
% fft for real-valued vectors
%
% use:
%   [b,f] = fftrl(a,t,mode)
%
% input:
%   a - input data
%   t - time vector
%   mode: 1: forward, -1:inverse
%

n = size(a);
a = reshape(a,[n(1) prod(n(2:end))]);

nt = length(t);
dt = t(2) - t(1);
nf = floor(nt/2) + 1;
tmax = t(end) - t(1);
f = 0:1/tmax:.5/dt;

switch mode
    case 1
        b = fft(a,[],1);
        b = b(1:nf,:);
        b = reshape(b,[nf n(2:end)]);
    case -1
        a = [a;conj(a(ceil(nt/2):-1:2,:))];
        b = ifft(a);
        b = real(b);
        b = reshape(b,[nt n(2:end)]);
    otherwise
        error('Unknown mode');
end

