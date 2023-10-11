function u = exprk3ds_cplx(u,F,g,A,m,tau,c2)

  t = 0;
  d = length(A{1});
  c3 = 2/3;

  eta1 = [7/4-3/sqrt(2)*1i,2/3-2/sqrt(3)*1i];
  alpha1 = [12/11+4*sqrt(2)/11*1i,3/4+sqrt(3)/4*1i];

  eta2 = 2^(d-2)*[-3+6*sqrt(2)*1i,-2/3+8/sqrt(3)*1i];
  alpha2 = [4/3+2*sqrt(2)/3*1i,6/7+3*sqrt(3)/7*1i];

  for mu = 1:d
    % First eq
    % For u2
    P1_1_u2{1}{mu} = phiquad(c2*tau*alpha1(1)*A{1}{mu},1);
    P2_1_u2{1}{mu} = phiquad(c2*tau*alpha2(1)*A{1}{mu},2);
    % For u3
    P1_1_u3{1}{mu} = phiquad(c3*tau*alpha1(1)*A{1}{mu},1);
    P2_1_u3{1}{mu} = phiquad(c3*tau*alpha2(1)*A{1}{mu},2);
    P1_2_u3{1}{mu} = phiquad(c3*tau*alpha1(2)*A{1}{mu},1);
    P2_2_u3{1}{mu} = phiquad(c3*tau*alpha2(2)*A{1}{mu},2);
    % For final approx 
    P1_1_u{1}{mu} = phiquad(tau*alpha1(1)*A{1}{mu},1);
    P2_1_u{1}{mu} = phiquad(tau*alpha2(1)*A{1}{mu},2);
    P1_2_u{1}{mu} = phiquad(tau*alpha1(2)*A{1}{mu},1);
    P2_2_u{1}{mu} = phiquad(tau*alpha2(2)*A{1}{mu},2);

    % Second eq
    % For u2
    P1_1_u2{2}{mu} = phiquad(c2*tau*alpha1(1)*A{2}{mu},1);
    P2_1_u2{2}{mu} = phiquad(c2*tau*alpha2(1)*A{2}{mu},2);
    % For u3
    P1_1_u3{2}{mu} = phiquad(c3*tau*alpha1(1)*A{2}{mu},1);
    P2_1_u3{2}{mu} = phiquad(c3*tau*alpha2(1)*A{2}{mu},2);
    P1_2_u3{2}{mu} = phiquad(c3*tau*alpha1(2)*A{2}{mu},1);
    P2_2_u3{2}{mu} = phiquad(c3*tau*alpha2(2)*A{2}{mu},2);
    % For final approx 
    P1_1_u{2}{mu} = phiquad(tau*alpha1(1)*A{2}{mu},1);
    P2_1_u{2}{mu} = phiquad(tau*alpha2(1)*A{2}{mu},2);
    P1_2_u{2}{mu} = phiquad(tau*alpha1(2)*A{2}{mu},1);
    P2_2_u{2}{mu} = phiquad(tau*alpha2(2)*A{2}{mu},2);
  end
  for jj = 1:m
    Fn{1} = F{1}(t,u{1},u{2});
    Fn{2} = F{2}(t,u{1},u{2});

    gn{1} = g{1}(t,u{1},u{2});
    gn{2} = g{2}(t,u{1},u{2});

    U2{1} = u{1} + c2*tau*real(eta1(1)*tucker(Fn{1},P1_1_u2{1}) + ...
                               eta2(1)*tucker(Fn{1},P2_1_u2{1}));
    U2{2} = u{2} + c2*tau*real(eta1(1)*tucker(Fn{2},P1_1_u2{2}) + ...
                               eta2(1)*tucker(Fn{2},P2_1_u2{2}));
    gn2diff{1} = g{1}(t+c2*tau,U2{1},U2{2})-gn{1};
    gn2diff{2} = g{2}(t+c2*tau,U2{1},U2{2})-gn{2};

    U3{1} = u{1} + c3*tau*real(eta1(1)*tucker(Fn{1},P1_1_u3{1}) + ...
                               eta2(1)*tucker(Fn{1},P2_1_u3{1})) + ...
                   4/9/c2*tau*real(eta1(2)*tucker(gn2diff{1},P1_2_u3{1}) + ...
                                   eta2(2)*tucker(gn2diff{1},P2_2_u3{1}));
    U3{2} = u{2} + c3*tau*real(eta1(1)*tucker(Fn{2},P1_1_u3{2}) + ...
                               eta2(1)*tucker(Fn{2},P2_1_u3{2})) + ...
                   4/9/c2*tau*real(eta1(2)*tucker(gn2diff{2},P1_2_u3{2}) + ...
                                   eta2(2)*tucker(gn2diff{2},P2_2_u3{2}));
    gn3diff{1} = g{1}(t+c3*tau,U3{1},U3{2})-gn{1};
    gn3diff{2} = g{2}(t+c3*tau,U3{1},U3{2})-gn{2};

    u{1} = u{1} + tau*real(eta1(1)*tucker(Fn{1},P1_1_u{1}) + ...
                           eta2(1)*tucker(Fn{1},P2_1_u{1})) + ...
                  3/2*tau*real(eta1(2)*tucker(gn3diff{1},P1_2_u{1}) + ...
                               eta2(2)*tucker(gn3diff{1},P2_2_u{1}));
    u{2} = u{2} + tau*real(eta1(1)*tucker(Fn{2},P1_1_u{2}) + ...
                           eta2(1)*tucker(Fn{2},P2_1_u{2})) + ...
                  3/2*tau*real(eta1(2)*tucker(gn3diff{2},P1_2_u{2}) + ...
                               eta2(2)*tucker(gn3diff{2},P2_2_u{2}));

    t = t + tau;
  end
end
