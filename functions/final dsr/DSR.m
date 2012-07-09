function CIG = DSR(data,t,xr,xs,z,v)

% form common offset image gather
DSR_step = opTEST([length(t),length(xr),length(xs)],t,xr,xs,z,v(1:length(z)));
image = DSR_step*vec(data); 

% shot_receiver -- midpoint_offset conversion
image = reshape(image,length(z),length(xr),length(xs));
MHImage = SR2MH(image,1); % conversion
[f,ky,kh,h] = kykh(t,xr,xs); % compute mh(yh) coordinate

% form common waveparameter gather via Radon
CIGmig = opCIGmig(f,ky,kh,h,z);
CIG = CIGmig*vec(MHImage);


