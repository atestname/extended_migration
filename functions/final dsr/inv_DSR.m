function data = inv_DSR(CIG,t,xr,xs,z,v)

% compute mh(yh) coordinate
[f,ky,kh,h] = kykh(t,xr,xs);

% form time midpoint offset image 
C = opCIGmig(f,ky,kh,h,z);
MHImage = C'*CIG;

% convert to shot receiver domain
MHImage = reshape(MHImage,length(z),length(ky),length(kh));
Image = SR2MH(MHImage,2);

% back propagate via inv DSR
DSR_step = opTEST([length(t),length(xr),length(xs)],t,xr,xs,z,v(1:length(z)));
data = DSR_step'*vec(Image);





