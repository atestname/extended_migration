function a = adjoint_test(K)

% perform the adjoint test for P

x_test = rand(size(K,2),1);%+rand(size(K,2),1)*i;
y_test = K*rand(size(K,2),1);%+rand(size(K,1),1)*i;

temp = (K*x_test);
left = y_test'*(K*x_test);
right = (K'*y_test)'*x_test;

a = norm(left-right)/max(size(K,2),size(K,1));


%adjoint_test(opDSR([10,10,10],t(1:10),xr(1:10),xs(1:10),z(1:10),v(1:10)));