% random jittering missing shots
n = ns;
p = .5;
I_jitter = jitter1d(n,p*n);
S_random = zeros(n,1); S_random(I_random) = 1;
S_jitter = zeros(n,1); S_jitter(I_jitter) = 1;
Js = opDiag(S_jitter);
Dt = opDirac(nt);
Dr = opDirac(nr);
RM = opKron(Js,Dr,Dt);



% random dithering in time
p = .5;
RM = opSimSourceRandTimeDither([nt,nr,ns],[p*nt*ns,1],ns);


% random weights of simu shots
p = .5;
nse = p*ns;
SS = rand(nse,ns);
Dt = opDirac(nt);
Dr = opDirac(nr);
RM = opKron(SS,Dr,Dt);

% simu data
simD = RM*D;

% sparsity representation
C = opCurvelet(nt, ns);

% operator
A = RM*C';
xest = spgl1(A,simD,0,1e-3,[],options);
f = C'*xest;
SNR = snr(D,f);

