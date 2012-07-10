function [b,f,k] = fktran(a,t,x,mode)
% f-k transform for real-values input.
%
% use:
%    [b,f,k] = fktran(a,t,x,mode)
%

nt = length(t);
nx = length(x);
dt = t(2)-t(1);
dx = x(2)-x(1);
xmax = x(end) - x(1);

k = -.5/dx:1/xmax:.5/dx;

switch mode
    case 1
        [b,f] = fftrl(a,t,1);
        b = ifft(b,[],2);
        b = circshift(b,[0 ceil(nx/2)-1]);
    case -1  
        b = circshift(a,[0 floor(nx/2)+1]);
        b = fft(b,[],2);
        b = fftrl(b,t,-1);
        b = real(b);
    otherwise
        error('unknown mode');
end