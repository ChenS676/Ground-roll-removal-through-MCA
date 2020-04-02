function  res2 = ctranspose(A)
A.adjoint = xor(A.adjoint,1);
res2 = A;