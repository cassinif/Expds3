function U = etd2rkds(U,g,A,m,tau)
  t = 0;
  d = length(A{1});
  for mu = 1:d
    [P2{1}{mu},P1{1}{mu}] = phiellm(tau*A{1}{mu},2);
    [P2{2}{mu},P1{2}{mu}] = phiellm(tau*A{2}{mu},2);
  end
  for jj = 1:m
    gn{1} = g{1}(t,U{1},U{2});
    gn{2} = g{2}(t,U{1},U{2});
    Fn{1} = kronsumv(U{1},A{1}) + gn{1};
    Fn{2} = kronsumv(U{2},A{2}) + gn{2};

    U2{1} = U{1} + tau*tucker(Fn{1},P1{1});
    U2{2} = U{2} + tau*tucker(Fn{2},P1{2});
    Up{1} = U2{1} + 2^(d-1)*tau*tucker(g{1}(t+tau,U2{1},U2{2})-gn{1},P2{1});
    U{2} = U2{2} + 2^(d-1)*tau*tucker(g{2}(t+tau,U2{1},U2{2})-gn{2},P2{2});
    U{1} = Up{1};

    t = t + tau;
  end
end
