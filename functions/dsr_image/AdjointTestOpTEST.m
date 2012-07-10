% adjoint test for the operator
dim = [10,10,10];
t = 0:.1:.9;
x = 0:10:90;
y = 0:10:90;
z = 0:10:90;
v = 1000*ones(size(z));
K = opTEST(dim,t,x,y,z,v);
adjoint_test(K)

