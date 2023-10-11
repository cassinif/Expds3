function u = exprk3ds_real(u,F,g,A,m,tau,c2)

  t = 0;
  c3 = 2/3;

  d = length(A{1});

  if d == 2
    eta1 = [-5/4,-4/3];
    alpha1{1} = [4/15*sqrt(10)+4/3,sqrt(33)/8+9/8];
    alpha1{2} = [-4/15*sqrt(10)+4/3,-sqrt(33)/8+9/8]; %beta1

    eta2 = [9,22/3];
    alpha2{1} = [2/9*sqrt(10)+16/9,3/22*sqrt(33)+3/2];
    alpha2{2} = [-2/9*sqrt(10)+16/9,-3/22*sqrt(33)+3/2]; %beta2
    for mu =1:d
      % First eq
      % For u2
      P1_1_u2{1}{mu} = phiquad(c2*tau*alpha1{mu}(1)*A{1}{mu},1);
      P2_1_u2{1}{mu} = phiquad(c2*tau*alpha2{mu}(1)*A{1}{mu},2);
      % For u3
      P1_1_u3{1}{mu} = phiquad(c3*tau*alpha1{mu}(1)*A{1}{mu},1);
      P2_1_u3{1}{mu} = phiquad(c3*tau*alpha2{mu}(1)*A{1}{mu},2);
      P1_2_u3{1}{mu} = phiquad(c3*tau*alpha1{mu}(2)*A{1}{mu},1);
      P2_2_u3{1}{mu} = phiquad(c3*tau*alpha2{mu}(2)*A{1}{mu},2);
      % For final approx 
      P1_1_u{1}{mu} = phiquad(tau*alpha1{mu}(1)*A{1}{mu},1);
      P2_1_u{1}{mu} = phiquad(tau*alpha2{mu}(1)*A{1}{mu},2);
      P1_2_u{1}{mu} = phiquad(tau*alpha1{mu}(2)*A{1}{mu},1);
      P2_2_u{1}{mu} = phiquad(tau*alpha2{mu}(2)*A{1}{mu},2);

      % Second eq
      % For u2
      P1_1_u2{2}{mu} = phiquad(c2*tau*alpha1{mu}(1)*A{2}{mu},1);
      P2_1_u2{2}{mu} = phiquad(c2*tau*alpha2{mu}(1)*A{2}{mu},2);
      % For u3
      P1_1_u3{2}{mu} = phiquad(c3*tau*alpha1{mu}(1)*A{2}{mu},1);
      P2_1_u3{2}{mu} = phiquad(c3*tau*alpha2{mu}(1)*A{2}{mu},2);
      P1_2_u3{2}{mu} = phiquad(c3*tau*alpha1{mu}(2)*A{2}{mu},1);
      P2_2_u3{2}{mu} = phiquad(c3*tau*alpha2{mu}(2)*A{2}{mu},2);
      % For final approx 
      P1_1_u{2}{mu} = phiquad(tau*alpha1{mu}(1)*A{2}{mu},1);
      P2_1_u{2}{mu} = phiquad(tau*alpha2{mu}(1)*A{2}{mu},2);
      P1_2_u{2}{mu} = phiquad(tau*alpha1{mu}(2)*A{2}{mu},1);
      P2_2_u{2}{mu} = phiquad(tau*alpha2{mu}(2)*A{2}{mu},2);
    end
    for jj = 1:m
      Fn{1} = F{1}(t,u{1},u{2});
      Fn{2} = F{2}(t,u{1},u{2});

      gn{1} = g{1}(t,u{1},u{2});
      gn{2} = g{2}(t,u{1},u{2});

      U2{1} = u{1} + c2*tau*(eta1(1)*tucker(Fn{1},P1_1_u2{1}) + ...
                             eta2(1)*tucker(Fn{1},P2_1_u2{1}));
      U2{2} = u{2} + c2*tau*(eta1(1)*tucker(Fn{2},P1_1_u2{2}) + ...
                             eta2(1)*tucker(Fn{2},P2_1_u2{2}));
      gn2diff{1} = g{1}(t+c2*tau,U2{1},U2{2})-gn{1};
      gn2diff{2} = g{2}(t+c2*tau,U2{1},U2{2})-gn{2};

      U3{1} = u{1} + c3*tau*(eta1(1)*tucker(Fn{1},P1_1_u3{1}) + ...
                             eta2(1)*tucker(Fn{1},P2_1_u3{1})) + ...
                     4/9/c2*tau*(eta1(2)*tucker(gn2diff{1},P1_2_u3{1}) + ...
                                 eta2(2)*tucker(gn2diff{1},P2_2_u3{1}));
      U3{2} = u{2} + c3*tau*(eta1(1)*tucker(Fn{2},P1_1_u3{2}) + ...
                             eta2(1)*tucker(Fn{2},P2_1_u3{2})) + ...
                     4/9/c2*tau*(eta1(2)*tucker(gn2diff{2},P1_2_u3{2}) + ...
                                 eta2(2)*tucker(gn2diff{2},P2_2_u3{2}));
      gn3diff{1} = g{1}(t+c3*tau,U3{1},U3{2})-gn{1};
      gn3diff{2} = g{2}(t+c3*tau,U3{1},U3{2})-gn{2};

      u{1} = u{1} + tau*(eta1(1)*tucker(Fn{1},P1_1_u{1}) + ...
                         eta2(1)*tucker(Fn{1},P2_1_u{1})) + ...
                    3/2*tau*(eta1(2)*tucker(gn3diff{1},P1_2_u{1}) + ...
                             eta2(2)*tucker(gn3diff{1},P2_2_u{1}));
      u{2} = u{2} + tau*(eta1(1)*tucker(Fn{2},P1_1_u{2}) + ...
                         eta2(1)*tucker(Fn{2},P2_1_u{2})) + ...
                    3/2*tau*(eta1(2)*tucker(gn3diff{2},P1_2_u{2}) + ...
                             eta2(2)*tucker(gn3diff{2},P2_2_u{2}));

      t = t + tau;
    end
  else
    eta1 = [2243/1350+440521/675/sqrt(2991111),19/27+151/27/sqrt(2391)];
    eta2 = 2^(d-3)*[-12544/675,-196/27];
    eta3 = [2243/1350-440521/675/sqrt(2991111),19/27-151/27/sqrt(2391)];
    alpha1 = [3*(5161+sqrt(2991111))/15869,3*(121+sqrt(2391))/490];
    alpha2 = [45/28,9/7];
    alpha3 = [3*(5161-sqrt(2991111))/15869,3*(121-sqrt(2391))/490];
    for mu =1:d
      % First eq
      % For u2
      P1_1_u2{1}{mu} = phiquad(c2*tau*alpha1(1)*A{1}{mu},1);
      P2_1_u2{1}{mu} = phiquad(c2*tau*alpha2(1)*A{1}{mu},2);
      P1b_1_u2{1}{mu} = phiquad(c2*tau*alpha3(1)*A{1}{mu},1);
      % For u3
      P1_1_u3{1}{mu} = phiquad(c3*tau*alpha1(1)*A{1}{mu},1);
      P2_1_u3{1}{mu} = phiquad(c3*tau*alpha2(1)*A{1}{mu},2);
      P1b_1_u3{1}{mu} = phiquad(c3*tau*alpha3(1)*A{1}{mu},1);
      P1_2_u3{1}{mu} = phiquad(c3*tau*alpha1(2)*A{1}{mu},1);
      P2_2_u3{1}{mu} = phiquad(c3*tau*alpha2(2)*A{1}{mu},2);
      P1b_2_u3{1}{mu} = phiquad(c3*tau*alpha3(2)*A{1}{mu},1);
      % For final approx 
      P1_1_u{1}{mu} = phiquad(tau*alpha1(1)*A{1}{mu},1);
      P2_1_u{1}{mu} = phiquad(tau*alpha2(1)*A{1}{mu},2);
      P1b_1_u{1}{mu} = phiquad(tau*alpha3(1)*A{1}{mu},1);
      P1_2_u{1}{mu} = phiquad(tau*alpha1(2)*A{1}{mu},1);
      P2_2_u{1}{mu} = phiquad(tau*alpha2(2)*A{1}{mu},2);
      P1b_2_u{1}{mu} = phiquad(tau*alpha3(2)*A{1}{mu},1);
      % Second eq
      % For u2
      P1_1_u2{2}{mu} = phiquad(c2*tau*alpha1(1)*A{2}{mu},1);
      P2_1_u2{2}{mu} = phiquad(c2*tau*alpha2(1)*A{2}{mu},2);
      P1b_1_u2{2}{mu} = phiquad(c2*tau*alpha3(1)*A{2}{mu},1);
      % For u3
      P1_1_u3{2}{mu} = phiquad(c3*tau*alpha1(1)*A{2}{mu},1);
      P2_1_u3{2}{mu} = phiquad(c3*tau*alpha2(1)*A{2}{mu},2);
      P1b_1_u3{2}{mu} = phiquad(c3*tau*alpha3(1)*A{2}{mu},1);
      P1_2_u3{2}{mu} = phiquad(c3*tau*alpha1(2)*A{2}{mu},1);
      P2_2_u3{2}{mu} = phiquad(c3*tau*alpha2(2)*A{2}{mu},2);
      P1b_2_u3{2}{mu} = phiquad(c3*tau*alpha3(2)*A{2}{mu},1);
      % For final approx 
      P1_1_u{2}{mu} = phiquad(tau*alpha1(1)*A{2}{mu},1);
      P2_1_u{2}{mu} = phiquad(tau*alpha2(1)*A{2}{mu},2);
      P1b_1_u{2}{mu} = phiquad(tau*alpha3(1)*A{2}{mu},1);
      P1_2_u{2}{mu} = phiquad(tau*alpha1(2)*A{2}{mu},1);
      P2_2_u{2}{mu} = phiquad(tau*alpha2(2)*A{2}{mu},2);
      P1b_2_u{2}{mu} = phiquad(tau*alpha3(2)*A{2}{mu},1);
    end
    for jj = 1:m
      Fn{1} = F{1}(t,u{1},u{2});
      Fn{2} = F{2}(t,u{1},u{2});

      gn{1} = g{1}(t,u{1},u{2});
      gn{2} = g{2}(t,u{1},u{2});

      U2{1} = u{1} + c2*tau*(eta1(1)*tucker(Fn{1},P1_1_u2{1}) + ...
                             eta2(1)*tucker(Fn{1},P2_1_u2{1}) + ...
                             eta3(1)*tucker(Fn{1},P1b_1_u2{1}));
      U2{2} = u{2} + c2*tau*(eta1(1)*tucker(Fn{2},P1_1_u2{2}) + ...
                             eta2(1)*tucker(Fn{2},P2_1_u2{2}) + ...
                             eta3(1)*tucker(Fn{2},P1b_1_u2{2}));
      gn2diff{1} = g{1}(t+c2*tau,U2{1},U2{2})-gn{1};
      gn2diff{2} = g{2}(t+c2*tau,U2{1},U2{2})-gn{2};

      U3{1} = u{1} + c3*tau*(eta1(1)*tucker(Fn{1},P1_1_u3{1}) + ...
                             eta2(1)*tucker(Fn{1},P2_1_u3{1}) + ...
                             eta3(1)*tucker(Fn{1},P1b_1_u3{1})) + ...
                     4/9/c2*tau*(eta1(2)*tucker(gn2diff{1},P1_2_u3{1}) + ...
                                 eta2(2)*tucker(gn2diff{1},P2_2_u3{1}) + ...
                                 eta3(2)*tucker(gn2diff{1},P1b_2_u3{1}));
      U3{2} = u{2} + c3*tau*(eta1(1)*tucker(Fn{2},P1_1_u3{2}) + ...
                             eta2(1)*tucker(Fn{2},P2_1_u3{2}) + ...
                             eta3(1)*tucker(Fn{2},P1b_1_u3{2})) + ...
                     4/9/c2*tau*(eta1(2)*tucker(gn2diff{2},P1_2_u3{2}) + ...
                                 eta2(2)*tucker(gn2diff{2},P2_2_u3{2}) + ...
                                 eta3(2)*tucker(gn2diff{2},P1b_2_u3{2}));
      gn3diff{1} = g{1}(t+c3*tau,U3{1},U3{2})-gn{1};
      gn3diff{2} = g{2}(t+c3*tau,U3{1},U3{2})-gn{2};

      u{1} = u{1} + tau*(eta1(1)*tucker(Fn{1},P1_1_u{1}) + ...
                         eta2(1)*tucker(Fn{1},P2_1_u{1}) + ...
                         eta3(1)*tucker(Fn{1},P1b_1_u{1})) + ...
                    3/2*tau*(eta1(2)*tucker(gn3diff{1},P1_2_u{1}) + ...
                             eta2(2)*tucker(gn3diff{1},P2_2_u{1}) + ...
                             eta3(2)*tucker(gn3diff{1},P1b_2_u{1}));
      u{2} = u{2} + tau*(eta1(1)*tucker(Fn{2},P1_1_u{2}) + ...
                         eta2(1)*tucker(Fn{2},P2_1_u{2}) + ...
                         eta3(1)*tucker(Fn{2},P1b_1_u{2})) + ...
                    3/2*tau*(eta1(2)*tucker(gn3diff{2},P1_2_u{2}) + ...
                             eta2(2)*tucker(gn3diff{2},P2_2_u{2}) + ...
                             eta3(2)*tucker(gn3diff{2},P1b_2_u{2}));

      t = t + tau;
    end
  end
end
