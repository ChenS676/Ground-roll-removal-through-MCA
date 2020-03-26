function  res1 = ctranspose(A)
A.adjoint = xor(A.adjoint,1);
res1 = A;