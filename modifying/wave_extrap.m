function v = wave_extrap(u,t,x,v,dz,dir)
% wavefield extrapolation by phase shift
%
% use:
%   v = wave_extrap(u,t,x,v,dz,dir)
%
% input:
%   u   - wavefield as matrix of size length(t) x length(x)
%   t   - time coordinates in seconds as column vector
%   x   - space coordinates in meters as row vector
%   v   - velocity in m/s (scalar)   
%   dz  - depth step in meters
%   dir - 1: forward in time, -1 backward in time
%
% output:
%   v   - extrapolated wavefield
%

% fk transform of the wavefield
[spec,f,k]= fktran(u,t,x,1);

% define grid and kz
[ff,kk]   = ndgrid(f,k);
kz        = 2*pi*sqrt((ff/v).^2-kk.^2);

% set sign of real part of kz
kz        = -sign(dir)*real(kz)+1i*abs(imag(kz));

% apply phase shift
spec      = exp(1i*abs(dz)*kz).*spec;

% inverse fk transform
v         = fktran(spec,t,k,-1);

% take real part of wavefield
v = real(v);