function v = DSR_step(u,ff,kkx,kky,dz,v)


F1 = opFFT3(size(ff));
F2 = opFFT2(size(ff));

% DSR operators
Px = 2*pi*sqrt((ff/v).^2-kkx.^2);
Py = 2*pi*sqrt((ff/v).^2-kky.^2);

Px = real(Px)+1i*abs(imag(Px));
Py = real(Py)+1i*abs(imag(Py));

E = exp(1i*abs(dz)*(Px + Py));


v = F2'*opDiag(vec(E))*F1*vec(u);

