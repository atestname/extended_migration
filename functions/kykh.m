function [f ky kh h] =  kykh(t,x,y)

n = length(x);
[xx yy] = ndgrid(x,y);
hh = .5*(yy - xx);
yy = .5*(yy + xx);
hmin = min(min(hh)); hmax = max(max(hh)); 
ymin = min(min(yy)); ymax = max(max(yy));
h = linspace(hmin,hmax,2*n);
y = linspace(ymin,ymax,n);
%h = floor(h);
y = floor(y);
[f, ky kh] = fkk(t,y,h);

