function  res=radon_transform1(s)
res.adjoint = 0;

% initialize parameters of radon operator
res.dt = s.dt;
res.nt = s.nt;
res.dx = s.dx;
res.nx = s.nx;
res.np = s.np;
res.pmax = s.pmax;

% initialize tau- and p series 
res.x = s.x;
res.p = s.p;
res.t = s.t;

res = class(res,'radon_transform1');