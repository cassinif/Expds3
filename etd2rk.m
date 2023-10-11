function u = etd2rk(u,Fv,gv,Kfun,m,tau,tol)
  t = 0;

  mkry1 = 1;
  mkry2 = 1;
  iom = 2;
  nn = length(u);
  for jj = 1:m
    Fn = Fv(t,u);
    gn = gv(t,u);
    [utmp,mkry1] = phipm_simul_iom(tau,Kfun,[zeros(nn,1),Fn],tol,mkry1,iom);
    U2 = u + utmp;
    [utmp2,mkry2] = phipm_simul_iom(tau,Kfun,[zeros(nn,2),1/tau*(gv(t+tau,U2)-gn)],tol,mkry2,iom);
    u = U2 + utmp2;
    t = t + tau;
  end
end
