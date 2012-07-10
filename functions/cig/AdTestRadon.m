% adjoint test for opRadon
t = 0:.1:1;
h = -100:10:100;
q = 1e-3*[-10:1:10];
power = 1;
op = opRadon(t,h,q,power);

adjoint_test(op)