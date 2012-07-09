clear;
close all;

% time coordinate in seconds as column vector
t = [0:.004:1]';

% spatial coordinate in meters as row vector
x = 0:10:1270;
z = 0:10:1270;

% generate 2D grid:
[tt,xx] = ndgrid(t,x);

% velocity
% background velocity [m/s]
v = 2000*ones(length(z));


figure;imagesc(v); colormap(gray);
% Source wavefield is impulsive source at t=0.1 and x=500,
source = (tt-.1).*exp(-1e-3*(xx-500).^2 - 1e3*(tt-.1).^2);
figure;imagesc(real(source));colormap(gray)


Ft = opDFT(length(t));
Ix = opDirac(length(x));
Fx = opDFT(length(x));
F = opKron(Fx,Ft);

D = F*vec(source);

data = D;
data = source; 

image = DSRImagingDFT2D(data,v(1:length(x)),z,x,t,2); 
%image = DSRImagingDFT3D(data,v(1:length(z)),z,x,model.freq,1); 

figure; imagesc(real(image))
