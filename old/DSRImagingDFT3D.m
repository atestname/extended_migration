function [image,uu] = DSRImagingDFT3D(u,c,z,x,tf,flag)
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
    F = opKron(Fs,Fr,Ft);
    U = F*vec(u);
    t = tf;
    fnyq = 1. / (2*(t(2)-t(1)));
    nf = size(u,1);
    %f = linspace(0.,fnyq,nf)';
    f = [0:df:fnyq-df -fnyq:df:-df];
end
% compute w 
% w = 2*pi*f;
w = f;
% compute kr(wavenumber for receiver)
dx = (x(2)-x(1));
fnyq = 1. / (2*(dx));
nf = size(u,2);
dff = 2*fnyq/nf;
fr = [0:dff:fnyq-dff -fnyq:dff:-dff];
% kr = fr;
% kx=[0:dkx:kxnyq-dkx -kxnyq:dkx:-dkx];
% kr = 2*pi*fr;
kr = sort(fr);
ks = kr;






% propagate wavefield in frequency-wavenumber domain using DSR with
% finite difference method 
% U = reshape(U,length(tf),length(x),length(x));
UU = zeros(length(tf)*length(x)*length(x),length(z));
UU0 = UU; %UU1 = UU; UU2 = UU; UU3 = UU; UU4 = UU;
DSR = zeros(length(tf),length(kr),length(ks));
%DSR = zeros(length(tf),length(kr));
dz = z(2)-z(1);
U0 = U; %U1 = U; U2 = U; U3 = U; U4 = U;
for iz = 1:length(z) % loop through all depth level
    % compute DSR for this depth level
    %keyboard;
    cc = c(iz);
    for iw = 1:length(w) 
        temp =  (w(iw)).^2;
        for ikr = 1:length(kr)
            SR2 = 2*pi*sqrt(temp./cc^2 - (kr(ikr))^2); % square root two
            SR2 = real(SR2)-1i*abs(imag(SR2));
            for iks = 1:length(ks)
                SR1 = 2*pi*sqrt(temp./cc^2 - (ks(iks))^2);
                SR1 = real(SR1)-1i*abs(imag(SR1));
                DSR(iw,ikr,iks) = -1j.*(SR1+SR2);  % Mr. DSR for this depth
%                 DSR1 = real(DSR)-1i*abs(imag(DSR));  % Mr. DSR for this depth
%                 DSR2 = real(DSR)+1i*abs(imag(DSR));  
%                 DSR3 = -real(DSR)-1i*abs(imag(DSR)); 
%                 DSR4 = -real(DSR)+1i*abs(imag(DSR));
            end
        end
    end
    % U = U + dz.*vec(DSR).*U;
    U0 = U0.*vec(exp(DSR.*dz));
%     U1 = U1.*vec(exp(DSR1.*dz));
%     U2 = U2.*vec(exp(DSR2.*dz));
%     U3 = U3.*vec(exp(DSR3.*dz));
%     U4 = U4.*vec(exp(DSR4.*dz));
    UU0(:,iz) = U0;
%     UU1(:,iz) = U1;
%     UU2(:,iz) = U2;
%     UU3(:,iz) = U3;
%     UU4(:,iz) = U4;
end % end of loop through all depth level

% transform U back to frequency space domain;
Iz = opDirac(length(z));
%Ft = opFFTsym_conv_mask(2*length(w)-1);
Ft = opDFT(length(w));
F_inv = opKron(Iz,Fs,Fr,Ft');
uu = F_inv*vec(UU0);
uu = reshape(uu,length(w),length(x),length(x),length(z));
%uut = sum(uu,1);
% form the image, sum all the frequencies getting zero time
for ix = 1:length(x) % loop through all horizontal place
    for iz = 1:length(z)
        image(iz,ix) = uu(1,ix,ix,iz);
    end
end
figure; imagesc(real(image))

% for ix = 1:length(x) % loop through all horizontal place
%     for iz = 1:length(z)
%         image(iz,ix) = sum(uu(1,:,ix,iz),2);
%     end
% end
% 
% figure; imagesc(real(image))



1+1;






    