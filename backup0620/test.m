% read data
[data,SuTraceHeaders,SuHeader]=ReadSu('data_ex3.su');

% source and receiver coordinates per trace
for k=1:length(SuTraceHeaders)
    xs(k)= SuTraceHeaders(k).SourceX;
    xr(k)= SuTraceHeaders(k).GroupX;
end

% note that the sample interval is given in miliseconds!
t = [0:SuHeader.ns-1]'*SuHeader.dt*1e-6;

% source and receiver coordinate vectors
xs = unique(xs);
xr = unique(xr);

% reshape data into cube
data = reshape(data,length(t),length(xr),length(xs));

% window and subsample
data(1:100,:,:) = 0;
data = data(:,1:10:end,:);
xr = xs;

% image using Tristan's code
z      = 0:10:1000;
v      = 2000;
image1 = DSR_mig(data,t,xr,xs,z,v);

figure;imagesc(image1);


% % image using Lina's code
% 
% image2 = my_mig(data,t,xr,z,v+0*z,2);
% 
% figure;imagesc(real(image2));