function y = dj_CEM_2(j, n_T, m, K, Delta, precision)

k_sum = 0;
for k = 0:j
    l_sum = 0;
    for l = 0:k
        n_sum = 0;
        % Evaluate only the iterations of the sum over n for which n >= k - 2*l
        for n = min(max(0,k-2*l),n_T):n_T
            switch precision
                case 'double'   % Double-precision arithmetic
                    contrib = ((((-K*Delta/2)^(2*l-k+2*n)) * gamma(j+m+2*l-k+2*n)) / (factorial(n) * ((m+K)^(j+m+2*l-k+2*n)) * gamma(2*l-k+n+1)));
                    if ~isnan(contrib)
                        n_sum = n_sum + contrib;
                    end
                case 'symbolic' % Symbolic arithmetic
                    n_sum = n_sum + ((((-K*Delta/2)^sym(2*l-k+2*n)) * gamma(sym(j+m+2*l-k+2*n))) / (factorial(sym(n)) * ((m+K)^sym(j+m+2*l-k+2*n)) * gamma(sym(2*l-k+n+1))));
                otherwise
                    error('Incorrect precision option')
            end
        end
        l_sum = l_sum + ( nchoosek(k,l) * n_sum );
    end
    k_sum = k_sum + ( (nchoosek(j,k) * ((Delta/2)^k)) * l_sum );
end

y = k_sum;
