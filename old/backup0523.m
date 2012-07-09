% a version not using kron and spot
% function [image] = DSRImaging(u,c,z,x,tf,flag)
% This dfunction is used to perform the DSR migration for lateral
% homogeneous media, basic formulation can be found at Equation 13 
% on http://sepwww.stanford.edu/sep/prof/bei/sg/paper_html/node11.html#eqn:3p9

% input
%   
%   - flag = 1, data is given in frequency domain(default, if not specified 
%               by user)
%            - u(f,r,s)
%            is the 3D data volumn of a 2D seismic survy observed at the
%            surface z = 0 in frequency domain, first dimension is frequency,
%            second is receiver, third is source 
%            - tf is the frequency coordinate 
% 
%   - flag = 2, data is given in time domain
%            - u(t,r,s) 
%            is the 3D data volumn of a 2D seismic survy observed at the
%            surface z = 0 in time domain, the first dimension is time, the
%            second is receiver, the third is source
%            - tf is the time coordinate
%
%   c(z)     -is the velocity vector as a function of z
%   z        -vertical coordinate--depth--(equally spaced) 
%   x        -horizontal coordinate--survey line--(equally spaced)
%             x is both source and receiver grid
%   
% output
%   image(z,x) -is the migrated image, a 2D data volum describing the
%              subsurface 


% initialize image
image = zeros(length(z),length(x));

if not(exist('flag','var'))
    flag = 1;
end

if flag == 1 %data is given in frequency domain
    U = u;
    for p = 2:3
        U = fft(U,[],p);
    end
    f = tf; 
else % data is given in time domain
    U = fftn(u);
    t = tf;
    fnyq = 1. / (2*(t(2)-t(1)));
    nf = size(U,1);
    %f = linspace(0.,fnyq,nf)';
    f = [linspace(0,fnyq,ceil(nf/2)),linspace(-fnyq,0,floor(nf/2))];
end
% compute w 
w = 2*pi*f;
% compute kr(wavenumber for receiver)
df = 1/(x(2)-x(1));
fnyq = 1. / (2*(df));
nf = size(U,2);
fr = [linspace(0,fnyq,ceil(nf/2)),linspace(-fnyq,0,floor(nf/2))];
kr = 2*pi*fr;
% compute ks(wavenumber for source)
% fnyq = 1. / (2*(x(2)-x(1)));
% nf = size(U,3);
% f = linspace(0.,fnyq,nf)';
% ks = 2*pi*f;
ks = kr; % as receiver and source are on the same grid



% propagate wavefield in frequency-wavenumber domain using DSR with
% finite difference method 
UU = zeros(size(U,1),size(U,2),size(U,3),length(z));
DSR = zeros(length(w),length(ks),length(kr));
dz = z(2)-z(1);
dz = 1;
for iz = 1:length(z) % loop through all depth level
    % compute DSR for this depth level
    for iw = 1:length(w) 
        for ikr = 1:length(kr)
            for iks = 1:length(ks)
                temp = w(iw)^2/c(iz)^2;
                SR1 = sqrt(temp - (ks(iks))^2); % square root one
                SR2 = sqrt(temp - (kr(ikr))^2); % square root two
                DSR(iw,ikr,iks) = -1j*(SR1+SR2);  % Mr. DSR for this depth
                %U(iw,ikr,iks) = (1+dz*DSR)*U(iw,ikr,iks); % propagate both source and receiver wavefield with DSR
            end
        end
    end
    U = U + dz.*DSR.*U;
    UU(:,:,:,iz) = U;
end % end of loop through all depth level
% transform U back to time space domain
uu = UU;
for p = 1:3
    uu = ifft(uu,[],p);
end
% if flag == 1 % data was in frequency domain
%     for p = 1:3
%         uu = ifftshift(uu,p);
%     end
% else % data was in time domain
     uu = ifftshift(uu);
% end
% form the image

for ix = 1:length(x) % loop through all horizontal place
    for iz = 1:length(z)
        image(iz,ix) = uu(iw,ix,ix,iz);
        
    end
end






    