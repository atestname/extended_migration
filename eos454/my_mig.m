function image = my_mig(data,t,xr,xs,z,v)
%
% migration
%
% use:
%   image = my_mig(data,t,xr,xs,z,v)
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
    source = (tt-.1).*exp(-1e-3*(xx-xs(is)).^2 - 1e3*(tt-.1).^2);
    receiver = data(:,:,is);
    % loop over depth levels
    % initialize shoti and reci
    shoti = source;
    reci = receiver;
    for iz = 1:length(z)
        % use my_step to advance one depthlevel
        shoti = my_step(shoti,t,xr,v,(iz-1)*dz,1);
        reci = my_step(reci,t,xr,v,(iz-1)*dz,-1);
        % update image 
        image(iz,:) = image(iz,:)+sum(shoti.*reci,1);
    end
    % end loop over depth levels
   
end
% end loop over shots



