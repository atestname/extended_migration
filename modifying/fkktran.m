function [b,f,kx,ky] = fkktran(a,t,x,y,mode)
% f-k-k transform for real-values input.
%
% use:
%    [b,f,kx,ky] = fktran(a,t,x,y,mode)
%

nt = length(t);
nx = length(x);
ny = length(y);
dt = t(2)-t(1);
dx = x(2)-x(1);
dy = y(2)-y(1);
xmax = x(end) - x(1);
ymax = y(end) - y(1);

kx = -.5/dx:1/xmax:.5/dx;
ky = -.5/dy:1/ymax:.5/dy;

switch mode
    case 1
        [b,f] = fftrl(a,t,1);
        b = ifft(b,[],2);
        b = ifft(b,[],3);
        b = circshift(b,[0 ceil(nx/2)-1 ceil(ny/2)-1]);
    case -1  
        b = circshift(a,[0 floor(nx/2)+1 floor(ny/2)+1]);
        b = fft(b,[],3);
        b = fft(b,[],2);
        b = fftrl(b,t,-1);
        b = real(b);
    otherwise
        error('unknown mode');
end