function image = mig_extrap(data,t,xr,xs,z,v)
%
% shot-receiver migration for constant velocity by wavefield extrapolation.
%
% use:
%   image = mig_extrap(data,t,xr,xs,v)
%
% input:
%   data  - data cube of size length(t) x length(xr) x length(xs) 
%   t     - time coordinate in seconds as column vector
%   xr    - receiver coordinate in meters as row vector
%   xs    - source coordinate in meters as row vector
%   z     - depth coordinate in meters as column vector
%   v     - velocity in m/s (scalar)
%
% output:
%   image - image as matrix of size length(z) x length(xr)

% initialize image
image = zeros(length(z),length(xr));

% depth step 
dz = z(2) - z(1);

% t-x grid
[tt,xx] = ndgrid(t,xr);

% loop over shots
for is = 1:length(xs)
    % construct source and receiver wavefields
    source   = (tt-.1).*exp(-1e-3*(xx-xs(is)).^2 - 1e3*(tt-.1).^2);
    
    % loop over depth levels
    for iz = 1:length(z)
        % use my_step to advance one depthlevel
        sourcei      = wave_extrap(source,t,xr,v,iz*dz,1);
        receiveri    = wave_extrap(data(:,:,is),t,xr,v,iz*dz,-1);
        
        % update image
        image(iz,:) = image(iz,:) + sum(sourcei.*receiveri,1);
    end
    % end loop over depth levels
end
% end loop over shots
