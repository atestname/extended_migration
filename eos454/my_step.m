function v = my_step(u,t,x,v,dz,dir)
% wavefield extrapolation by phase shift
%
% use:
%   v = my_step(u,t,x,v,dz,dir)
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

[spec,f,k]= fktran(u,t,x);
[ff,kk]   = ndgrid(f,k);

kz        = 2*pi*sqrt((ff/v).^2-kk.^2);

kz        = -sign(dir)*real(kz)+1i*abs(imag(kz));

spec      = exp(1i*abs(dz)*kz).*spec;

v         = ifktran(spec,f,k);
v         = v(1:length(t),1:length(x));