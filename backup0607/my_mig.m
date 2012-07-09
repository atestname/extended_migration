function image = my_mig(data,tf,x,z,v,flag)
%
% migration
%
% use:
%   image = my_mig(data,t,xr,xs,z,v)
%
% input:
%   u    -
%       flag = 1, data is given in frequency domain(default, if not specified 
%                   by user)
%                - u(f,r,s)
%                is the 3D data volumn of a 2D seismic survy observed at the
%                surface z = 0 in frequency domain, first dimension is frequency,
%                second is receiver, third is source 
%                - tf is the frequency coordinate 
%     
%       flag = 2, data is given in time domain
%                - u(t,r,s) 
%                is the 3D data volumn of a 2D seismic survy observed at the
%                surface z = 0 in time domain, the first dimension is time, the
%                second is receiver, the third is source
%                - tf is the time coordinate
%   x     - receiver and coordinate in meters as row vector
%   z     - depth coordinate in meters as column vector
%   v     - velocity in m/s (scalar)
%
% output:
%   image - image as matrix of size length(z) x length(xr)

% initialize image
image = zeros(length(z),length(x));

% depth step 
dz = z(2) - z(1);


for iz = 1:length(z)
    ppi = my_step(data,tf,x,v(iz),(iz-1)*dz,1);
    pp = reshape(ppi,length(tf),length(x),length(x),flag);
    % update image
    for ix = 1:length(x)
        image(iz,ix) = image(iz,ix)+pp(1,ix,ix);
    end
% ppp = permute(pp,[2,3,1]);
% figure(2);imagesc(real(pp(:,:,33)));colormap(gray);figure(3);imagesc(real(image))
% colorbar;zlim([1000 5000]);
end




