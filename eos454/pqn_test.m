
% Number of time samples
nt = 1024;
% Number of sources
ns = 178;
% Number of receivers
nr = 178;

% Time sampling interval
dt = 0.004;

% Read data
D = ReadSuFast('GulfOfSuez178.su');
D = reshape(D,nt,nr,ns);

% Select small subset
D = D(1:1024,30,1:175);

% Define new data sizes
[nt,nr,ns] = size(D);

% Vectorize D
D = D(:);

% % Display
% figure
% imagesc(reshape(D,nt,ns)); colormap(gray); colorbar;
% title('Original data (receiver gather)');
% xlabel('Shot number'); ylabel('Time sample')

C = opCurvelet(nt, ns);

% Use the spgl1 and pqnl1 as inverse and display the result. 
% fid = fopen('log.txt', 'w'); 
options = spgSetParms('optTol', 1e-4, 'iterations', 200);%, 'fid', fid); 
p_list = [.3];
for i = 1:length(p_list)
    p = p_list(i);
    D_RM = opSimSourceRandTimeDither([nt,nr,ns],[p*nt*ns,1],ns);
    simD = D_RM*D;
    A = D_RM*C';
    [xestspg, ~,~, infospg] = spgl1(A,simD,0,1e-3,[],options);
    infospg.snr = snr(D,C'*xestspg(:,i));
    Destspg = C'*xestspg(:,i);
    opts.iterations = 50; opts.optTol = 1e-4;opts.verbosity = 1;
    [xestpqn, ~,~, infopqn] = pqnl1(A,simD,0,1e-3,[],opts);
    %[xestpqn,infopqn] = bpdn_pqn(A,simD,1e-3,0,opts);
    infopqn.snr = snr(D,C'*xestpqn(:,i));
    Destpqn = C'*xestpqn(:,i);
   
    

%     f(:,i) = C'*xest(:,i);
%     SNR(i) = snr(D,f(:,i)); 
end

fprintf('info spg')
infospg
fprintf('info pqn')
infopqn


     

timespg = infospg.timeTotal;
timepqn = infopqn.timeTotal;
    

figure;plot(1,timespg,1,timepqn);legend('spg','pqn')

%
figure('Name','Solution paths')
plot(infospg.xNorm1,infospg.rNorm2,infopqn.xNorm1,infopqn.rNorm2)
hold on;
scatter(infospg.xNorm1,infospg.rNorm2);
scatter(infopqn.xNorm1,infopqn.rNorm2);
%scatter(infoamp.xNorm1,infoamp.rNorm2);
hold off
xlabel('one-norm model')
ylabel('two-norm residual')
title('Solutions paths')
legend('SPGl1','PQNl1')
axis tight


figure; 
subplot(1,2,1);imagesc(reshape(Destspg,nt,ns)); colormap(gray);
title(strcat(['SPGL1,p = .3, SNR=' num2str(infospg.snr) 'dB']))
subplot(1,2,2);imagesc(reshape(Destpqn,nt,ns)); colormap(gray);
title(strcat(['PQNl1,p = .3, SNR=' num2str(infopqn.snr) 'dB']))
