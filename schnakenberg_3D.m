clear all
close all

d = 3;

n = 80*ones(1,d);
ld = 0*ones(1,d);
rd = 1*ones(1,d);
T = 0.6;

deltau = 0.05;
deltav = 1;
alphau = 0.01;
alphav = 0.01;
rho = 100;
au = 0.1305;
av = 0.7695;

tol_kry = 1e-6;

mrange_etd2rkds = 800:200:1400;
mrange_exprk3ds_cplx = 200:100:500;
mrange_exprk3ds_real = 200:100:500;
mrange_etd2rk = 800:200:1400;
mrange_exprk3 = 200:100:500;

c2 = 1/3; % for exprk3

for mu = 1:d
  x{mu} = linspace(ld(mu),rd(mu),n(mu));
  h(mu) = (rd(mu)-ld(mu))/(n(mu)-1);
  D2{mu} = spdiags(ones(n(mu),1)*([1,-2,1]/(h(mu)^2)),-1:1,n(mu),n(mu));
  D2{mu}(1,1:2) = [-2,2]/(h(mu)^2);
  D2{mu}(n(mu),(n(mu)-1):n(mu)) = [2,-2]/(h(mu)^2);
  D1{mu} = spdiags(ones(n(mu),1)*([-1,0,1]/(2*h(mu))),-1:1,n(mu),n(mu));
  D1{mu}(1,1:2) = [0,0];
  D1{mu}(n(mu),(n(mu)-1):n(mu)) = [0,0];
  A_sp{1}{mu} = -alphau*D1{mu} + deltau*D2{mu};
  A_sp{2}{mu} = -alphav*D2{mu} + deltav*D2{mu};
  A{1}{mu} = full(A_sp{1}{mu});
  A{2}{mu} = full(A_sp{2}{mu});
end
[X{1:d}] = ndgrid(x{1:d});

g{1} = @(t,u,v) rho*(au-u+(u.*u).*v);
g{2} = @(t,u,v) rho*(av-(u.*u).*v);

F{1} = @(t,u,v) kronsumv(u,A{1}) + g{1}(t,u,v);
F{2} = @(t,u,v) kronsumv(v,A{2}) + g{2}(t,u,v);

K{1} = kronsum(A_sp{1});
K{2} = kronsum(A_sp{2});

pn = prod(n);
Fvec = @(t,u) [K{1}*u(1:pn)+g{1}(t,u(1:pn),u(pn+1:2*pn));...
              K{2}*u(pn+1:2*pn)+g{2}(t,u(1:pn),u(pn+1:2*pn))];

Kfun = @(uvec) [K{1}*uvec(1:pn);K{2}*uvec(pn+1:2*pn)];

gvec = @(t,uvec) [g{1}(t,uvec(1:pn),uvec(pn+1:2*pn));g{2}(t,uvec(1:pn),uvec(pn+1:2*pn))];

U0{1} = au+av+1e-5*exp(-100*((X{1}-1/3).^2+(X{2}-1/2).^2+(X{3}-1/3).^2));
U0{2} = (av/(au+av)^2)*ones(n);
u0 = [U0{1}(:);U0{2}(:)];

%nstepsref = 2000;
%tau = T/nstepsref;
disp('Computing/loading reference solution...')
%Uref = exprk3ds_real(U0,F,g,A,nstepsref,tau,c2);
load('ref_SCHNAK_3D.mat','Uref')
disp('Reference solution computed/loaded!')
uref = [Uref{1}(:),Uref{2}(:)];
normref = norm(uref,inf);

counter = 0;
fprintf('ETD2RK\n')
for m = mrange_etd2rk
  fprintf('m=%i\n',m)
  counter = counter + 1;
  tau = T/m;
  tic
  Uetd2rk = etd2rk(u0,Fvec,gvec,Kfun,m,tau,tol_kry);
  cpu_etd2rk(counter) = toc;
  err_etd2rk(counter) = norm(Uetd2rk(:)-uref(:),inf)/normref;
end
counter = 0;

fprintf('EXPRK3\n')
for m = mrange_exprk3
  fprintf('m=%i\n',m)
  counter = counter + 1;
  tau = T/m;
  tic
  Urk3 = exprk3(u0,Fvec,gvec,Kfun,m,tau,tol_kry,c2);
  cpu_exprk3(counter) = toc;
  err_exprk3(counter) = norm(Urk3(:)-uref(:),inf)/normref;
end

counter = 0;
fprintf('ETD2RKds\n')
for m = mrange_etd2rkds
  fprintf('m=%i\n',m)
  counter = counter + 1;
  tau = T/m;
  tic
  Uetd2rkds = etd2rkds(U0,F,g,A,m,tau);
  cpu_etd2rkds(counter) = toc;
  err_etd2rkds(counter)= norm([Uetd2rkds{1}(:);Uetd2rkds{2}(:)]-uref(:),inf)/normref;
end

counter = 0;
fprintf('EXPRK3ds (cplx)\n')
for m = mrange_exprk3ds_cplx
  fprintf('m=%i\n',m)
  counter = counter + 1;
  tau = T/m;
  tic
  Uexprk3ds_cplx = exprk3ds_cplx(U0,F,g,A,m,tau,c2);
  cpu_exprk3ds_cplx(counter) = toc;
  err_exprk3ds_cplx(counter)= norm([Uexprk3ds_cplx{1}(:);Uexprk3ds_cplx{2}(:)]-uref(:),inf)/normref;
end

counter = 0;
fprintf('EXPRK3ds (real)\n')
for m = mrange_exprk3ds_real
  fprintf('m=%i\n',m)
  counter = counter + 1;
  tau = T/m;
  tic
  Uexprk3ds_real = exprk3ds_real(U0,F,g,A,m,tau,c2);
  cpu_exprk3ds_real(counter) = toc;
  err_exprk3ds_real(counter)= norm([Uexprk3ds_real{1}(:);Uexprk3ds_real{2}(:)]-uref(:),inf)/normref;
end

figure;
loglog(mrange_etd2rkds,err_etd2rkds,'og',...
       mrange_exprk3ds_cplx,err_exprk3ds_cplx,'dm',...
       mrange_exprk3ds_real,err_exprk3ds_real,'+c',...
       mrange_etd2rk, err_etd2rk,'xr',...
       mrange_exprk3,err_exprk3,'sb',...
       mrange_etd2rkds,err_etd2rkds(end)*(mrange_etd2rkds/mrange_etd2rkds(end)).^(-2),'--k',...
       mrange_exprk3,err_exprk3(end)*(mrange_exprk3/mrange_exprk3(end)).^(-3),'-k')

legend('ETD2RKds','EXPRK3ds (cplx)','EXPRK3ds (real)','ETD2RK','EXPRK3')
title('Error decay')
xlabel('Number of time steps')
ylabel('Error in the infinity norm')
drawnow

figure;
loglog(cpu_etd2rkds,err_etd2rkds,'og',...
       cpu_exprk3ds_cplx,err_exprk3ds_cplx,'dm',...
       cpu_exprk3ds_real,err_exprk3ds_real,'+c',...
       cpu_etd2rk,err_etd2rk,'xr',...
       cpu_exprk3,err_exprk3,'sb')
legend('ETD2RKds','EXPRK3ds (cplx)','EXPRK3ds (real)','ETD2RK','EXPRK3')
title('CPU diagram')
xlabel('Wall-clock time (s)')
ylabel('Error in the infinity norm')
drawnow
