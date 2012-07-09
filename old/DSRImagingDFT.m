function [image] = DSRImagingDFT(u,c,z,x,tf,flag)
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

% transform data to w-ks-kr domain
% operators
Ft = opDFT(size(u,1));
Fr = opDFT(size(u,2));
Fs = opDFT(size(u,3));
if flag == 1 %data is given in frequency domain
    It = opDirac(length(tf));
    F = opKron(Fs',Fr',It);
    U = F*vec(u);
    f = tf; 
else % data is given in time domain
    F = opKron(Fs',Fr',Ft);
    U = F*vec(u);
    t = tf;
    fnyq = 1. / (2*(t(2)-t(1)));
    nf = size(u,1);
    %f = linspace(0.,fnyq,nf)';
    f = [0:df:fnyq-df -fnyq:df:-df];
end
% compute w 
w = 2*pi*f;
% compute kr(wavenumber for receiver)
df = 1/(x(2)-x(1));
fnyq = 1. / (2*(df));
nf = size(u,2);
dff = 2*fnyq/nf;
fr = [0:dff:fnyq -fnyq:dff:-dff];
% kx=[0:dkx:kxnyq-dkx -kxnyq:dkx:-dkx];
kr = 2*pi*fr;
kr = fr;
% compute ks(wavenumber for source)
% fnyq = 1. / (2*(x(2)-x(1)));
% nf = size(U,3);
% f = linspace(0.,fnyq,nf)';
% ks = 2*pi*f;
ks = kr; % as receiver and source are on the same grid



% propagate wavefield in frequency-wavenumber domain using DSR with
% finite difference method 
% U = reshape(U,length(tf),length(x),length(x));
UU = zeros(length(tf)*length(x)*length(x),length(z));
DSR = zeros(length(tf),length(ks),length(kr));
dz = z(2)-z(1);
for iz = 1:length(z) % loop through all depth level
    % compute DSR for this depth level
    for iw = 1:length(w) 
        for ikr = 1:length(kr)
            for iks = 1:length(ks)
                SR1 = sqrt(1./c(iz)^2 - (ks(iks)./w(iw))^2); % square root one
                SR2 = sqrt(1./c(iz)^2 - (kr(ikr)./w(iw))^2); % square root two
                DSR(iw,ikr,iks) = -1j.*w(iw).*(SR1+SR2);  % Mr. DSR for this depth
            end
        end 
        
    end
    % U = U + dz.*vec(DSR).*U;
    DSR1 = +real(DSR)-1i*abs(imag(DSR));
    DSR2 = -real(DSR)-1i*abs(imag(DSR));
    DSR3 = -real(DSR)+1i*abs(imag(DSR));
    DSR4 = +real(DSR)+1i*abs(imag(DSR));
    U = U.*vec(exp(DSR2.*dz));
    UU(:,iz) = U;
end % end of loop through all depth level

% transform U back to frequency space domain;
Iz = opDirac(length(z));
Ft = opFFTsym_conv_mask(2*length(w)-1);
F_inv = opKron(Iz,Fs,Fr,Ft');
uu = F_inv*vec(UU);
uu = reshape(uu,2*length(w)-1,length(x),length(x),length(z));
%uut = sum(uu,1);
% form the image, sum all the frequencies getting zero time
for ix = 1:length(x) % loop through all horizontal place
    for iz = 1:length(z)
        image(iz,ix) = uu(1,ix,ix,iz);
    end
end


figure; imagesc(real(image))




    