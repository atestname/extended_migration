function [ coefs,f,k ] = my_fktran( data,t,x )
%my_fktran I wanted to really understand the Fourier transform, and so
%wrote my own. No offense to CREWES. This one keeps the positive and
%negative frequencies because it's easier to understand that way, and for
%ease of inversion.

% Sampling period
dt = abs(t(2)-t(1));
% Sampling frequency
fs = 1/dt;
% Nyquist frequency
fnyq = fs/2;
% Number of samples
Nt = length(t);
% Spatial equivalents
dx = abs(x(2)-x(1));
ks = 1/dx;
knyq = ks/2;
Nx = length(x);

f = fftshift((0:Nt-1)'*fs/Nt);
% The greater than cutoffs below were chosen to be consistent with
% fftshift, which puts the Nyquist frequency at the far negative end of the
% frequency interval.
f(f>=fnyq) = f(f>=fnyq)-fs;
k = fftshift((0:Nx-1)*ks/Nx);
k(k>=knyq) = k(k>=knyq)-ks;

coefs = fftshift(fft2(data));

end

