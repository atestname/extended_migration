function data = inv_DSR(CIG,t,xr,xs,z,v)

% compute mh(yh) coordinate
[ky,kh,h] = kykh(t,xr,xs);

% form time midpoint offset image 
MHImage = opCIGmig(t,ky,kh,h,z)'*CIG;

% convert to shot receiver domain
Image = SR2MH(MHImage,2);

% back propagate via inv DSR
data = DSR_mig_inv_spot_2(Image,t,xr,xs,z,v);




