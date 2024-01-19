function varargout = phiellm(A,ell)

% Compute phi_j(A) for 0 <= mu <= ell
% The output is given in the order phi_{ell},phi_{ell-1},...,phi_0

d = 7; % degree of pade

s = max(0,ceil(1+log2(norm(A,inf))));

A = A/2^s;
ID = eye(size(A));

i = d;
Ncoeff = sum(cumprod([1,d-(0:i-1)])./cumprod([1,2*d+ell-(0:i-1)]).*(-1).^(0:i)./(factorial((0:i)).*factorial(ell+i-(0:i))));
Dcoeff = ((-1)^(i))*prod([(d-i+1):1:d])/(prod([(2*d+ell-i+1):1:(2*d+ell)])*factorial(i));
varargout{1} = Ncoeff;
D = Dcoeff;

for i = (d-1):-1:0
  Ncoeff = sum(cumprod([1,d-(0:i-1)])./cumprod([1,2*d+ell-(0:i-1)]).*(-1).^(0:i)./(factorial((0:i)).*factorial(ell+i-(0:i))));
  varargout{1} = A*varargout{1} + Ncoeff*ID ;
  Dcoeff = ((-1)^(i))*prod([(d-i+1):1:d])/(prod([(2*d+ell-i+1):1:(2*d+ell)])*factorial(i));
  D = A*D + Dcoeff*ID ;
end

varargout{1} = full(D\varargout{1});

if ell > 0
  fac = factorial(ell-1);
  for i = (ell-1):-1:1
    varargout{ell-i+1} = A*varargout{ell-i}+ID/fac;
    fac = fac/i;
  end
  varargout{ell+1} = A*varargout{ell} + ID;
end

% Squaring
for i = 1:s
  for j = ell:-1:1
    varargout{ell-j+1} = varargout{ell+1}*varargout{ell-j+1} + varargout{ell-j+1};
    fac = 1;
    for k = (j-1):-1:1
      varargout{ell-j+1} = varargout{ell-j+1} + varargout{ell-k+1}/fac;
      fac = fac*(j-k+1);
    end
    varargout{ell-j+1} = varargout{ell-j+1}/(2^j);
  end
  varargout{ell+1} = varargout{ell+1}*varargout{ell+1};
end
