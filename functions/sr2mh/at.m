% perform the adjoint test for P

x_test = rand(10,10,10);
y_test = rand(10,10,20);

left = vec(y_test)'*vec(SR2MH(x_test,1));
right = vec(SR2MH(y_test,2))'*vec(x_test);

a = norm(left-right)/max(size(K,2),size(K,1));
