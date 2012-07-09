function CIG = DSR(data,t,xr,xs,z,v)

% form common offset image gather
[image] = DSR_mig_spot_2(data,t,xr,xs,z,v(1:length(z)));

% shot_receiver -- midpoint_offset conversion
MHImage = SR2MH(image,1); % conversion
[ky,kh,h] = kykh(t,xr,xs); % compute mh(yh) coordinate

% form common waveparameter gather via Radon
CIG = opCIGmig(t,ky,kh,h,z)*MHImage;


