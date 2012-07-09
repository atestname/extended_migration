function output = lpradon(input,t,h,q,power,mode)
% linear and parabolic Radon transform and adjoint.
%
% use:
%   output = lpradon(input,t,h,q,power,mode)
%
% input:
%   input - input matrix or vector (length(t) x length(h) for forward, length(t) x length(q) for adjoint)
%   t     - time vector in seconds 
%   h     - offset vecror in meters
%   q     - radon parameter
%   power - 1: linear radon, 2: parabolic radon
%   mode  - 1: forward, -1: adjoint
%
% Wrapper for code provided by M. Sacchi, by Tristan van Leeuwen

fmin = 0.01;
fmax = 100;

if mode == 1
    %input  = reshape(input,length(t),length(h));
    output = inverse_radon_freq(input,t(2)-t(1),h,q,power,fmin,fmax,.1,'adj');
    %output = output(:);
else
    %input  = reshape(input,length(t),length(q));
    output = forward_radon_freq(input,t(2)-t(1),h,q,power,fmin,fmax);
    %output = output(:);
end

end

function [m] = inverse_radon_freq(d,dt,h,q,N,flow,fhigh,mu,sol)
%INVERSE_RADON_FREQ: Inverse linear or parabolic Radon transform. 
%                    Frequency domain alg.
%
%  [m] = inverse_radon_freq(d,dt,h,q,N,flow,fhigh,mu,sol)
% 
%  IN   d:     seismic traces   
%       dt:    sampling in sec
%       h(nh): offset or position of traces in meters
%       q(nq): ray parameters  if N=1
%              residual moveout at far offset if N=2
%       N:     N=1 Linear tau-p  
%              N=2 Parabolic tau-p
%       flow:  freq.  where the inversion starts in HZ (> 0Hz)
%       fhigh: freq.  where the inversion ends in HZ (< Nyquist) 
%       mu:    regularization parameter 
%       sol:   sol='ls' least-squares solution
%              sol='adj' adjoint
%
%  OUT  m:     the linear or parabolic tau-p panel
%
%  Reference: Hampson, D., 1986, Inverse velocity stacking for multiple elimination
%             Journal of the CSEG, vol 22, no 1., 44-55.
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%



 [nt,nh] = size(d);
 nq = max(size(q));

 %if N==2; h=h/max(abs(h));end;  
 nfft = 2*(2^nextpow2(nt));

 D = fft(d,nfft,1);
 M = zeros(nfft,nq);
 i = sqrt(-1);

 ilow  = floor(flow*dt*nfft)+1; 
 if ilow < 2; ilow=2; end;
 ihigh = floor(fhigh*dt*nfft)+1;
 if ihigh > floor(nfft/2)+1; ihigh=floor(nfft/2)+1; end

 Q = eye(nq)*nh;

 ilow = max(ilow,2);

 for ifreq=ilow:ihigh

  f = 2.*pi*(ifreq-1)/nfft/dt;
  L = exp(i*f*(h.^N)'*q); 
  y = D(ifreq,:)';
  xa = L'*y;
  A = L'*L + mu*Q;

 if isequal(sol,'ls');  
  x  =  A\xa; 
 end; 

 if isequal(sol,'adj');  
  x  =  xa; 
 end

 M(ifreq,:) = x';
 M(nfft+2-ifreq,:) = conj(x)';

end

 M(nfft/2+1,:) = zeros(1,nq);
 %m = real(ifft(M,[],1));
 m = (ifft(M,[],1));
 m = m(1:nt,:);
return;
end

function [d]=forward_radon_freq(m,dt,h,p,N,flow,fhigh)
%FORWARD_RADON_FREQ: Forward linear and parabolic Radon transform.
%                    Freq. domain algorithm
%
%  [d] = forward_radon_freq(m,dt,h,p,N,flow,fhigh);
%
%  IN   m:     the Radon panel, a matrix m(nt,np)
%       dt:    sampling in sec
%       h(nh): offset or position of traces in mts
%       p(np): ray parameter  to retrieve if N=1
%              curvature of the parabola if N=2
%       N:     N=1 linear tau-p
%              N=2 parabolic tau-p
%       flow:  min freq. in Hz
%       fhigh: max freq. in Hz
%
%  OUT  d:     data
%
%  Reference: Hampson, D., 1986, Inverse velocity stacking for multiple elimination,
%             Journal of the CSEG, vol 22, no 1., 44-55.
%
%
%  Copyright (C) 2008, Signal Analysis and Imaging Group
%  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
%  Author: M.D.Sacchi
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%


[nt,nq] = size(m);
nh = length(h);

 %if N==2;  
 %  h=(h/max(abs(h))).^2; 
 %end;   

 nfft = 2*(2^nextpow2(nt));

 M = fft(m,nfft,1);
 D = zeros(nfft,nh);
 i = sqrt(-1);

 ilow  = floor(flow*dt*nfft)+1; 
 if ilow < 1; ilow=1; end;
 ihigh = floor(fhigh*dt*nfft)+1;
 if ihigh > floor(nfft/2)+1; ihigh = floor(nfft/2)+1; end
 
 for ifreq=ilow:ihigh
  f = 2.*pi*(ifreq-1)/nfft/dt;
  L = exp(i*f*((h.^N)'*p  ));
  x = M(ifreq,:)';
  y = L * x; 
  D(ifreq,:) = y';
  D(nfft+2-ifreq,:) = conj(y)';
 end

 D(nfft/2+1,:) = zeros(1,nh);
 %d = real(ifft(D,[],1));
 d = (ifft(D,[],1));
 d = d(1:nt,:);

return;



end