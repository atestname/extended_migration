function [H,dH] = Helm2D(f,m,o,d,n,nb,beta,l,order)
% 2D Helmholtz operator, discretized using 13-point stencil with PML
% [Hustedt et. al. GJI 2004 (157) pp. 1269--1296, eq. (A4)]
%
% The PML functions are of the form: xi = 1 + i*(beta/w)*[0:1/(nb-1):1].^l
%
% Tristan van Leeuwen, 2012
% tleeuwen@eos.ubc.ca
%
% use:
%   [H,dH] = Helm2D(f,m,o,d,n,nb,{beta},{l})
%
% input:
%   f       - frequency [1/s]
%   m       - gridded squared slowness [s^2/m^2]
%   {o,d,n} - grid definition: z = o(1) + [0:n(1)-1]*d(1), etc.
%   nb      - number of points to use as absorbing boundary in each direction
%   {beta,l}- PML parameters, default beta = 100, l = 2.
%   order   - 2 or 4. default = 4.
%
% output:
%   H       - sparse Helmholtz matrix of size n(1)*n(2) x n(1)*n(2)
%  dH       - Jacobian matrix dH/dm

% input checking
if nargin < 7
    beta = 100;
end
if nargin < 8    
    l    = 2;
end
if nargin < 9
    order = 4;
end

%% total size and angular frequency
N  = prod(n);
w  = 2*pi*f;

%% PML, layer is 2 gridpoints bigger than domain to facilitate averaging
p1 = [1 - (beta/w)*1i*linspace(1,0,nb(1)+1)'.^l;  ones(n(1)-2*nb(1),1);  1 - (beta/w)*1i*linspace(0,1,nb(1)+1)'.^l];
p2 = [1 - (beta/w)*1i*linspace(1,0,nb(2)+1).^l    ones(1,n(2)-2*nb(2))   1 - (beta/w)*1i*linspace(0,1,nb(2)+1).^l];

p1 = repmat(p1 ,1     ,n(2)+2); p1 = p1(:);
p2 = repmat(p2 ,n(1)+2,1   );   p2 = p2(:);

%% averaging matrices, again 2 gridpoints bigger

% (+/-) 0.5 gridpoints
A1 = kron(speye(n(2)+2),spdiags(ones(n(1)+2,1)*[.5 .5],[0 1],n(1)+2,n(1)+2)); 
A2 = kron(spdiags(ones(n(2)+2,1)*[.5 .5],[0 1],n(2)+2,n(2)+2),speye(n(1)+2)); 

% (+/-) 1.5 gridpoints
B1 = kron(speye(n(2)+2),spdiags(ones(n(1)+2,1)*[.5 .5],[1 2],n(1)+2,n(1)+2));
B2 = kron(spdiags(ones(n(2)+2,1)*[.5 .5],[1 2],n(2)+2,n(2)+2),speye(n(1)+2));

%% Restriction matrix to get rid of ghostpoints
R  = kron(spdiags(ones(n(2)+1,1),1,n(2),n(2)+2),spdiags(ones(n(1)+1,1),1,n(1),n(1)+2));

%% 13 point stencil
%                    a1mmmm
%                    a1mm  
%                    a1m  
%       a2mmm-a2m-a2m-ac-a2p-a2pp-a2ppp
%                    a1p
%                    a1pp    
%                    a1ppp
%
% recover a second order scheme by using: 
if order==2
    c1 = 0;c2 = 1;c3 = 0;
else
    c1 = (1/24)^2; c2 = (9/8)^2; c3 = (9/8)*(1/24);
end

a1p   =  (c2./(A1*p1) + c3./(A1'*p1) + c3./(B1*p1) +  0./(B1'*p1))./(d(1).^2*p1);
a1m   =  (c3./(A1*p1) + c2./(A1'*p1) +  0./(B1*p1) + c3./(B1'*p1))./(d(1).^2*p1);

a1pp  = -(c3./(A1*p1) +  0./(A1'*p1) + c3./(B1*p1) +  0./(B1'*p1))./(d(1).^2*p1);
a1mm  = -( 0./(A1*p1) + c3./(A1'*p1) +  0./(B1*p1) + c3./(B1'*p1))./(d(1).^2*p1);

a1ppp =  ( 0./(A1*p1) +  0./(A1'*p1) + c1./(B1*p1) +  0./(B1'*p1))./(d(1).^2*p1);
a1mmm =  ( 0./(A1*p1) +  0./(A1'*p1) +  0./(B1*p1) + c1./(B1'*p1))./(d(1).^2*p1);

a2p   =  (c2./(A2*p2) + c3./(A2'*p2) + c3./(B2*p2) +  0./(B2'*p2))./(d(2).^2*p2);
a2m   =  (c3./(A2*p2) + c2./(A2'*p2) +  0./(B2*p2) + c3./(B2'*p2))./(d(2).^2*p2);

a2pp  = -(c3./(A2*p2) +  0./(A2'*p2) + c3./(B2*p2) +  0./(B2'*p2))./(d(2).^2*p2);
a2mm  = -( 0./(A2*p2) + c3./(A2'*p2) +  0./(B2*p2) + c3./(B2'*p2))./(d(2).^2*p2);

a2ppp =  ( 0./(A2*p2) +  0./(A2'*p2) + c1./(B2*p2) +  0./(B2'*p2))./(d(2).^2*p2);
a2mmm =  ( 0./(A2*p2) +  0./(A2'*p2) +  0./(B2*p2) + c1./(B2'*p2))./(d(2).^2*p2);

ac    = -a1p-a1m-a2p-a2m-a1pp-a1mm-a2pp-a2mm-a1ppp-a1mmm-a2ppp-a2mmm;

%% assemble sparse matrix
H  = w.^2*spdiags(m,0,N,N)...
    + spdiags(R*[a2mmm a2mm a2m a1mmm a1mm a1m ac a1p a1pp a1ppp a2p a2pp a2ppp],[-3*n(1) -2*n(1) -n(1) -3 -2 -1 0 1 2 3 n(1) 2*n(1) 3*n(1)],N,N);

%% Jacobian matrix
dH = w.^2*speye(N);





