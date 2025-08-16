function [RFM_j_T, CEM_j_T, CEM_n_T] = calculate_number_of_terms(m, K, Delta, target_trunc_error, show_log)

%-------------------------------------------------------------------------------
%  NUMBER OF ITERATIONS OF RFM - Using epsilon_T expression
%-------------------------------------------------------------------------------

if show_log
    fprintf('\nNUMBER OF ITERATIONS OF RFM - Using epsilon_T expression\n')
end
j_T = -1;
j_sum = 0;
current_trunc_error = Inf;
while current_trunc_error > target_trunc_error
    j_T = j_T+1;
    contrib = ((K^j_T * dj(j_T, m, K, Delta, 'double'))/factorial(j_T));
    j_sum = j_sum + contrib;
    current_trunc_error = 1 - ((m^m)/gamma(m))*j_sum;
    if show_log
        fprintf('  j_T = %d | trunc. error = %g\n', j_T, current_trunc_error)
    end
end

RFM_j_T = j_T;

%-------------------------------------------------------------------------------
%  NUMBER OF ITERATIONS OF CEM - Evaluating CDF numerically
%-------------------------------------------------------------------------------

if show_log
    fprintf('\nNUMBER OF ITERATIONS OF CEM - Evaluating CDF numerically\n')
end
j_T = 0;
n_T = -1;
EbN0 = 1;
avg_gamma = 1;
current_trunc_error = Inf;
trunc_error_vector_CEM_I = [];
while current_trunc_error > target_trunc_error
    
    n_T = n_T+1;
    
    % Evaluate CDF numerically -------------------------------------------------
    u_max = 50*avg_gamma;
    nof_points = 2e3;
    u = linspace(0, u_max, nof_points);
    total_sum = 0;
    for j = 0:j_T
        k_sum = 0;
        for k = 0:j
            l_sum = 0;
            for l = 0:k
                n_sum = 0;
                for n = 0:n_T
                    contrib = ((-u*K*Delta/2).^(2*l-k+2*n)) / (factorial(n)*gamma(2*l-k+n+1));
                    % contrib(isnan(contrib)) = 0;
                    n_sum = n_sum + contrib;
                end
                f_exact = (u.^(j+m-1)).*exp(-(m+K)*u).*besseli(2*l-k, -u*K*Delta); % This is the exact expression with n --> infinity
                f = (u.^(j+m-1)).*exp(-(m+K)*u).*n_sum;
                if abs(f(end)) > max(abs(f))*0.005
                    figure
                    plot(u, f_exact, 'b-', u, f, 'r--')
                    error('The range of values of u is not long enough')
                end
                % l_sum = l_sum + (nchoosek(k,l) * sum(((f(1:end-1)+f(2:end))/2)*(u(2)-u(1))));
                contrib = trapz(u,f);
                if ~isnan(contrib)
                    l_sum = l_sum + (nchoosek(k,l) * contrib);
                end
            end
            k_sum = k_sum + (nchoosek(j,k)*((Delta/2)^k)*l_sum);
        end
        total_sum = total_sum + (((K^j)/factorial(j))*k_sum);
    end
    cdf_y_num = ((m^m)/gamma(m)) * total_sum;
    % --------------------------------------------------------------------------
    
    current_trunc_error = 1 - cdf_y_num;
    trunc_error_vector_CEM_I = [trunc_error_vector_CEM_I current_trunc_error];
    if show_log
        % fprintf('  j_T = %d | n_T = %d | trunc. error = %g [%g]\n', j_T, n_T, current_trunc_error, trunc_error_vector_RFM(j_T+1))
        fprintf('  j_T = %d | n_T = %d | trunc. error = %g\n', j_T, n_T, current_trunc_error)
    end
    if current_trunc_error < 0 || ...
        ((numel(trunc_error_vector_CEM_I) > 1) && (trunc_error_vector_CEM_I(end-1) - trunc_error_vector_CEM_I(end) <= eps))
        if show_log
            fprintf('  === increasing j_T = %d ---> %d\n', j_T, j_T+1)
            % pause
        end
        if j_T == 0
            % This was the first iteration for j_T = 0. We take the last
            % value of n_T obtained here for j_T = 0 as the final value for
            % n_T, because if we wait until the last iteration of j_T then
            % n_T will always be zero and that value doesn't work. See the
            % comment provided a few lines below.
            CEM_n_T = n_T;
            j_T = RFM_j_T; % Jump to this value to speed up the search since we know that j_T for CEM will be greater than for RFM
        end
        j_T = j_T + 1;
        n_T = -1;
        current_trunc_error = Inf;
        trunc_error_vector_CEM_I = [];
    end
end

CEM_j_T = j_T;
% CEM_n_T = n_T; % The last n_T will always be n_T=0 and the plots
% show that this value does not provide an accurate version of PDF/CDF
