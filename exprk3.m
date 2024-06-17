function u = exprk3(u,gv,Kfun,m,tau,tol,c2)
  t = 0;
  c3 = 2/3;
  mkry1 = 1;
  mkry2 = 1;
  mkry3 = 1;
  iom = 2;
  nn = length(u);
  for jj = 1:m
    gn = gv(t,u);
    Fn = Kfun(u) + gn;

    [utmp,mkry1] = phipm_simul_iom(c2*tau,Kfun,[zeros(nn,1),Fn],tol,mkry1,iom);
    U2 = u + utmp;
    [utmp,mkry2] = phipm_simul_iom(c3*tau,Kfun,[zeros(nn,1),Fn,4/9/c2/c3^2/tau*(gv(t+c2*tau,U2)-gn)],tol,mkry2,iom);
    U3 = u + utmp;
    [utmp,mkry3] = phipm_simul_iom(tau,Kfun,[zeros(nn,1),Fn,3/2/tau*(gv(t+c3*tau,U3)-gn)],tol,mkry3,iom);
    u = u + utmp;

    t = t + tau;
  end
end
