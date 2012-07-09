% play with spg_mmv or spg_bpdn to see how this joint sparsity solver works
%% problem setting
m = 120; n = 512; k = 20; g = 10;% m rows, n cols, k nonzeros, 10 unknown vectors
X = zeros(n,g); p = randperm(n);
for i = 1:g
    X(:,i) = zeros(n,1); X(p(1:k),i) = i.*sign(randn(k,1));
end

A  = randn(m,n); [Q,R] = qr(A',0);  A = Q';
    
B  = A*X + 0.005 * randn(m,g);


%%
opts = spgSetParms('optTol',1e-4);
sigma = 1e-3;
[X_hat,R,G,INFO] = spg_mmv(A,B,sigma,opts);
