function res = mtimes(A,x)

if A.adjoint == 0                      % A*x
    x=reshape(x,A.nt,A.np);
    ress=invfwd_tx_sstackn_linear(x,A.dt,A.p,A.x);  
    res=reshape(ress,A.nt*A.nx,1);
else                                            % A.t*x
    x=reshape(x,A.nt,A.nx);
    ress=fwd_tx_sstackn_linear(x,A.dt,A.p,A.x);     
    res=reshape(ress,A.nt*A.np,1);
end