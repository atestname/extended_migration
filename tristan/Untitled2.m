model.o = [0 0];
model.d = [10 10];
model.n = [101 101];
model.nb = [30 30 0];
model.t0 = -.08;

t    = 0:.004:1; nt = length(t);
freq = [0:1/t(end):.5/(t(2)-t(1)), -.5/(t(2)-t(1)):1/t(end):-1/t(end)];
If   = [2:length(freq)];

model.freq = freq(If);             nfreq = length(model.freq);
model.zsrc = 10;
model.xsrc = linspace(0,1000,101);  nsrc  = length(model.xsrc);
model.zrec = 10;
model.xrec = linspace(0,1000,101); nrec  = length(model.xrec);
model.f0   = 10;
model.t0   = -.08;

Q = speye(nsrc); 

z = 0:10:1000;
x = 0:10:1000;
[zz,xx] = ndgrid(z,x);

% background velocity [m/s]
v0 = 2000 + 0.*zz;
dv = 0*xx;
dv(zz >= 520) = 600;
dv(zz >= 600) = 0;
figure; imagesc(dv);
v = v0+dv;
m = 1e6./v(:).^2;

[~,J] = F(m,Q,model);
Dobs = J*m;
save D Dobs
% A    = opKron(opRestriction(length(freq),If),opDirac(nrec),opDirac(nsrc));
% F =  A*opFFT1([nt,nrec,nsrc]);

F = opFFT1([nt,nrec,nsrc]);
Dt   = F'*Dobs;
xr = model.xrec; xs = model.xsrc; 

cube = permute(reshape(gather(Dt),nrec,nsrc,nt),[3 1 2]); 

a = test1(cube,t,xr,xs,z,v0(1:length(z)));
aa = DSR_mig(cube,t,xr,xs,z,v0(1:length(z)));

figure;imagesc(model.xrec,t,cube(:,:,26));colormap(seiscol(0,1));
