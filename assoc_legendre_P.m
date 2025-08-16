function y = assoc_legendre_P(mu, v, z, precision)
% Associated Legendre function of the first kind: P_v^mu(z)
%
% See Gradshteyn 8.7-8.8 Associated Legendre Functions (page 958)
% 
% Validation source:  with https://keisan.casio.com/exec/system/1180573406

% In the particular case considered in this problem, the argument z should
% always be real and greater than 1 given that m+K > K*Delta (m>0, K>0 and
% Delta lies within the interval [0,1]). The only way z could be a complex
% number is if we have that m+K < K*Delta, which leads to the square root
% of a negative number and therefore to a complex number, but this should
% never be the case. We double check that this never happens - just in case.
% if ~isreal(z)
%     error('Legendre function: The argument of the Legendre function should be real, but it is not!')
% end
% if z <= 1 || real(z) <= 1 || abs(z) <= 1
%     error('Legendre function: Input z is expected to be real and strictly greater than 1')
% end
% % Check input range
% if abs(1-z) >= 2
%     error('Legendre function: Input z is out of the expected range')
% end

switch precision
    case 'double'   % Double-precision arithmetic
        % This is the default option. Nothing special to be done here.
    case 'symbolic' % Symbolic arithmetic
        mu = sym(mu);
        v  = sym(v);
        z  = sym(z);
    otherwise
        error('Incorrect precision option')
end

if mu > 0 && mod(mu,1) == 0 % mu = 1, 2, 3, ... (non-zero positive integer)
    % The definition of the Legendre function of the first kind involves the
    % Gauss hypergeometric function 2F1, which loses its meaning when the third
    % input parameter is c = 1-mu = 0, -1, -2, -3, ... (non-interger positive),
    % in other words, when mu = 1, 2, 3, ... (non-zero positive integer). In
    % this case we need to use the alternative form provided in Gradshteyn
    % equation 9.101.1 (page 1005), which is an analytic continuation.
    % y = (((z+1)/(z-1))^(mu/2)) * ((pochhammer(-v,mu)*pochhammer(v+1,mu))/factorial(mu)) * (((1-z)/2)^mu) * hypergeom([mu-v, mu+v+1], mu+1, (1-z)/2);
    y = ((((1-z)/2) * sqrt((z+1)/(z-1)))^mu) * ((pochhammer(-v,mu)*pochhammer(v+1,mu))/factorial(mu)) * hypergeom([mu-v, mu+v+1], 1+mu, (1-z)/2);
else
    % Gradshteyn (8.702)
    y = (((z+1)/(z-1))^(mu/2)) * (1/gamma(1-mu)) * hypergeom([-v, v+1], 1-mu, (1-z)/2);
end
