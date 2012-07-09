function spec = fft3(u,dim)
if size(size(u))~=3
    u = reshape(u,dim);
end
spec = fft(u,[],1);
spec = fft(spec,[],2);
spec = fft(spec,[],3);
