function [f,ky,kh] = yh_convert(t,xr,xs)

% this function is to convert wavefield from shot offset domain to midpoint
% offset domain, according to the formulation 24 25 given at:
% http://sepwww.stanford.edu/sep/prof/bei/sg/paper_html/node12.html#SECTION00134000000000000000
% ky = .5*(ks-kh); => ky = ks + kg;
% kg = .5*(ks+kh); => kh = kg - ks;

% input
%     E : wavefield(f,ks,kg)
%     xr: receiver coordinate
%     xs: source coordinate

% output
%     F : wavefield(f,ks,kh)
%     ky: midpoint wavenumber
%     kh: half offset wavenumber
 
% compute ks(source wavenumber) kg(receiver wavenumber)
[f kg ks] = fkk(t,xr,xs);

% f-kg-ks grid
[kkg kks] = ndgrid(kg,ks);

% compute kh ky
% ky = .5*(ks-kh); => ky = ks + kg;
% kg = .5*(ks+kh); => kh = kg - ks;
kky = kks + kkg;
kkh = kkg - kks;

ky = unique(abs(kky));
kh = unique(abs(kkh));


% rotate the matrix 45' counter clockwise

