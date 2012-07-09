%Sequential-source data reconstruction (acquistion with randomly jittered missing shots)
% Number of time samples
nt = 1024;
% Number of sources
ns = 178;
% Number of receivers
nr = 178;

% Read data
D = ReadSuFast('../data/GulfOfSuez178.su');
D = reshape(D,nt,nr,ns);

% Select small subset
D = D(1:256,30,1:50);

% Define new data sizes
[nt,nr,ns] = size(D);

% Vectorize D
D = D(:);

% Display
figure
imagesc(reshape(D,nt,ns)); colormap(gray); colorbar;
title('Original data (receiver gather)');
xlabel('Shot number'); ylabel('Time sample')
