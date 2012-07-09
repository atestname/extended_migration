function [mdata] = zeromig(data1_zo,datazoTHeaders,datazoHeader,velocity)
% extract the time, source, receiver coordinates
% t0 = datazoHeader.dtOrig;
t0 = 0;
dt = datazoHeader.dt;
nt = size(data1_zo,1);
t = t0+[0:1:nt-1].*dt; % time coordinate
t = t./10^6; % time coordinate in second

source = [datazoTHeaders.SourceX];
receiver = [datazoTHeaders.GroupX];

% show zero offset profile in depth 
zz = .5.*t.*velocity;
figure;imagesc(source,zz,data1_zo); colormap(gray);
title('zero offset profile for data2_zo')
% fk transform for data
[spec,f,kx] = fktran(data1_zo,t,source,0,0,0,0);
% filter 5-40 frequency
j = 1;
for i = 1:length(f)
    if 5<=f(i) && f(i)<=40
        ff(j) = f(i);
        spec_filtered(j,:) = spec(i,:);
        j = j+1;
    end
end
% zero offset migration
z = 0:5:1000;
parms = .5*velocity*ones(length(z),length(source));
mdata = fx_zero_mig(spec_filtered,ff,parms,5,5,40);
figure; imagesc(source,z,mdata); title('migration data2_zo'); colormap(gray);

