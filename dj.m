function y = dj(j, m, K, Delta, precision)

outer_sum = 0;
for k = 0:j
    inner_sum = 0;
    for l = 0:k
        switch precision
            case 'double'   % Double-precision arithmetic
                  % Expression in the paper by Jiayi Zhang
                  % inner_sum = inner_sum + ( nchoosek(k,l) * gamma(j+m+2*l-k) * ((((m+K)^2)-((K*Delta)^2))^(-(j+m)/2)) * ...
                  %             exp(pi*(2*l-k)*1i/2) * assoc_legendre_P(k-2*l, j+m-1, (m+K)/sqrt(((m+K)^2)-((K*Delta)^2)), 'double') );
                  % Corrected expression
                  inner_sum = inner_sum + ( nchoosek(k,l) * gamma(j+m+2*l-k) * ((((m+K)^2)-((K*Delta)^2))^(-(j+m)/2)) * ...
                              ((-1)^(2*l-k)) * assoc_legendre_P(k-2*l, j+m-1, (m+K)/sqrt(((m+K)^2)-((K*Delta)^2)), 'double') );
            case 'symbolic' % Symbolic arithmetic
                  % Expression in the paper by Jiayi Zhang
                  % inner_sum = inner_sum + ( nchoosek(sym(k),sym(l)) * gamma(sym(j+m+2*l-k)) * ((((m+K)^2)-((K*Delta)^2))^(-(sym(j+m))/2)) * ...
                  %             exp(pi*(sym(2*l-k))*1i/2) * assoc_legendre_P(sym(k-2*l), sym(j+m-1), (m+K)/sqrt(((m+K)^2)-((K*Delta)^2)), 'symbolic') );
                  % Corrected expression
                  inner_sum = inner_sum + ( nchoosek(sym(k),sym(l)) * gamma(sym(j+m+2*l-k)) * ((((m+K)^2)-((K*Delta)^2))^(-(sym(j+m))/2)) * ...
                              (sym(-1)^sym(2*l-k)) * assoc_legendre_P(sym(k-2*l), sym(j+m-1), (m+K)/sqrt(((m+K)^2)-((K*Delta)^2)), 'symbolic') );
            otherwise
                error('Incorrect precision option')
        end
    end
    outer_sum = outer_sum + ( (nchoosek(j,k) * ((Delta/2)^k)) * inner_sum );
end

y = outer_sum;
