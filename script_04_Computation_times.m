
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE COMPUTATION TIMES IN THIS PC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if version ~=  '9.2.0.556344 (R2017a)'
    error('This program should be run in Matlab R2017a (version 9.2.0.556344)')
end

clear
clc
warning off

% Load configuration parameters and other data
load('Number_of_terms_and_iterations.mat', 'm', 'K', 'Delta', ...
                                           'target_trunc_error_values', ...
                                           'RFM_j_T_values', ...
                                           'CEM_j_T_values', ...
                                           'CEM_n_T_values')

% Other configuration parameters
precision = 'double';
nof_repetitions_RFM = 10;
nof_repetitions_CEM = 100;

% Vectors for results
comp_time_RFM_values      = zeros(nof_repetitions_RFM, numel(target_trunc_error_values));
comp_time_CEM_1_values    = zeros(nof_repetitions_CEM, numel(target_trunc_error_values));
comp_time_CEM_2_values    = zeros(nof_repetitions_CEM, numel(target_trunc_error_values));
comp_time_CEM_3_values    = zeros(nof_repetitions_CEM, numel(target_trunc_error_values));

% Calculate computation times - Do not use parallel computing (parfor)
for i = 1:numel(target_trunc_error_values)
    
    target_trunc_error = target_trunc_error_values(i);
    RFM_j_T = RFM_j_T_values(i);
    CEM_j_T = CEM_j_T_values(i);
    CEM_n_T = CEM_n_T_values(i);
    
    fprintf('\nCalculating results for eps_T = %.2e ( %d | %d | %d ) [%d/%d]: ', target_trunc_error, RFM_j_T, CEM_j_T, CEM_n_T, i, numel(target_trunc_error_values))
    
    % -------------------------------------------------------------------------
    % Evaluation of RFM
    % -------------------------------------------------------------------------
    fprintf('RFM...')
    for r = 1:nof_repetitions_RFM
        this_d_j = 0;
        tic
        for j = 0:RFM_j_T
            this_d_j = this_d_j + dj(j, m, K, Delta, precision);
        end
        comp_time_RFM_values(r, i) = toc;
    end
    
    % -------------------------------------------------------------------------
    % Evaluation of CEM-I
    % -------------------------------------------------------------------------
    fprintf('CEM1...')
    for r = 1:nof_repetitions_CEM
        this_d_j = 0;
        tic
        for j = 0:CEM_j_T
            this_d_j = this_d_j + dj_CEM_1(j, CEM_n_T, m, K, Delta, precision);
        end
        comp_time_CEM_1_values(r, i) = toc;
    end
    
    % -------------------------------------------------------------------------
    % Evaluation of CEM-II
    % -------------------------------------------------------------------------
    fprintf('CEM2...')
    for r = 1:nof_repetitions_CEM
        this_d_j = 0;
        tic
        for j = 0:CEM_j_T
            this_d_j = this_d_j + dj_CEM_2(j, CEM_n_T, m, K, Delta, precision);
        end
        comp_time_CEM_2_values(r, i) = toc;
    end
    
    % -------------------------------------------------------------------------
    % Evaluation of CEM-III
    % -------------------------------------------------------------------------
    fprintf('CEM3...')
    for r = 1:nof_repetitions_CEM
        this_d_j = 0;
        tic
        for j = 0:CEM_j_T
            this_d_j = this_d_j + dj_CEM_3(j, CEM_n_T, m, K, Delta, precision);
        end
        comp_time_CEM_3_values(r, i) = toc;
    end
    
    fprintf('OK')
    
end
fprintf('\n\n')

figure
loglog(target_trunc_error_values, mean(comp_time_RFM_values,1), 'k-', ...
       target_trunc_error_values, mean(comp_time_CEM_1_values,1), 'r-', ...
       target_trunc_error_values, mean(comp_time_CEM_2_values,1), 'g-', ...
       target_trunc_error_values, mean(comp_time_CEM_3_values,1), 'b-')
set(gca, 'XLim', [target_trunc_error_values(end) target_trunc_error_values(1)])
xlabel('Truncation error (epsilon\_T)')
ylabel('Computation time (sec)')
legend('RFM', 'CEM-I', 'CEM-II', 'CEM-III')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE FLOPS IN THIS PC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% IMPORTANT - Read the information in the following links:
% https://uk.mathworks.com/help/matlab/matlab_prog/measure-performance-of-your-program.html
% https://uk.mathworks.com/help/matlab/matlab_prog/profiling-for-improving-performance.html
% https://uk.mathworks.com/matlabcentral/answers/25962-processing-time-of-multiplicatoin-and-addition

fprintf('\nEstimating computation time for ADDITION ... ')
nof_operations = 1e9;
nof_repetitions = 5e2;
x_values = linspace(0, 50, nof_repetitions);
y_values = linspace(0, 50, nof_repetitions);
comp_times = zeros(1, nof_repetitions);
for i = 1:nof_repetitions
    x = x_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    y = y_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    result = 0; % Pre-allocate memory to store the result (avoids having to find free memory while timing)
    tic
    for j = 1:nof_operations
        % This is the operation being timed
        result = x + y;
        % ---------------------------------
    end
    comp_times(i) = toc/nof_operations;
end
comp_times_w_a = comp_times;
fprintf('OK\n')

fprintf('\nEstimating computation time for MULTIPLICATION ... ')
nof_operations = 1e9;
nof_repetitions = 5e2;
x_values = linspace(0, 50, nof_repetitions);
y_values = linspace(0, 50, nof_repetitions);
comp_times = zeros(1, nof_repetitions);
for i = 1:nof_repetitions
    x = x_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    y = y_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    result = 0; % Pre-allocate memory to store the result (avoids having to find free memory while timing)
    tic
    for j = 1:nof_operations
        % This is the operation being timed
        result = x * y;
        % ---------------------------------
    end
    comp_times(i) = toc/nof_operations;
end
comp_times_w_m = comp_times;
fprintf('OK\n')

fprintf('\nEstimating computation time for DIVISION ... ')
nof_operations = 1e9;
nof_repetitions = 5e2;
x_values = linspace(0, 50, nof_repetitions);
y_values = linspace(0, 50, nof_repetitions);
comp_times = zeros(1, nof_repetitions);
for i = 1:nof_repetitions
    x = x_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    y = y_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    result = 0; % Pre-allocate memory to store the result (avoids having to find free memory while timing)
    tic
    for j = 1:nof_operations
        % This is the operation being timed
        result = x / y;
        % ---------------------------------
    end
    comp_times(i) = toc/nof_operations;
end
comp_times_w_d = comp_times;
fprintf('OK\n')

fprintf('\nEstimating computation time for EXPONENTIATION ... ')
nof_operations = 1e7;
nof_repetitions = 5e2;
x_values = linspace(0, 100, nof_repetitions);
y_values = linspace(0, 100, nof_repetitions);
comp_times = zeros(1, nof_repetitions);
for i = 1:nof_repetitions
    x = x_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    y = y_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    result = 0; % Pre-allocate memory to store the result (avoids having to find free memory while timing)
    tic
    for j = 1:nof_operations
        % This is the operation being timed
        result = x ^ y;
        % ---------------------------------
    end
    comp_times(i) = toc/nof_operations;
end
comp_times_w_e = comp_times;
fprintf('OK\n')

fprintf('\nEstimating computation time for SQUARE ROOT ... ')
nof_operations = 1e8;
nof_repetitions = 5e2;
x_values = linspace(0, 100, nof_repetitions);
comp_times = zeros(1, nof_repetitions);
for i = 1:nof_repetitions
    x = x_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    result = 0; % Pre-allocate memory to store the result (avoids having to find free memory while timing)
    tic
    for j = 1:nof_operations
        % This is the operation being timed
        result = sqrt(x);
        % ---------------------------------
    end
    comp_times(i) = toc/nof_operations;
end
comp_times_w_r = comp_times;
fprintf('OK\n')

fprintf('\nEstimating computation time for FACTORIAL ... ')
nof_operations = 1e6;
nof_repetitions = 5e2;
x_values = round(linspace(0, 15, nof_repetitions));
comp_times = zeros(1, nof_repetitions);
for i = 1:nof_repetitions
    x = x_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    result = 0; % Pre-allocate memory to store the result (avoids having to find free memory while timing)
    tic
    for j = 1:nof_operations
        % This is the operation being timed
        result = factorial(x);
        % ---------------------------------
    end
    comp_times(i) = toc/nof_operations;
end
comp_times_w_f = comp_times;
fprintf('OK\n')

fprintf('\nEstimating computation time for BINOMIAL COEFFICIENT ... ')
nof_operations = 1e6;
nof_repetitions = 5e2;
x_values = round(linspace(0, 50, nof_repetitions));
y_values = round(linspace(0, 25, nof_repetitions));
comp_times = zeros(1, nof_repetitions);
for i = 1:nof_repetitions
    x = x_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    y = y_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    result = 0; % Pre-allocate memory to store the result (avoids having to find free memory while timing)
    tic
    for j = 1:nof_operations
        % This is the operation being timed
        result = nchoosek(x,y);
        % ---------------------------------
    end
    comp_times(i) = toc/nof_operations;
end
comp_times_w_b = comp_times;
fprintf('OK\n')

fprintf('\nEstimating computation time for GAMMA FUNCTION ... ')
nof_operations = 1e7;
nof_repetitions = 5e2;
x_values = linspace(0, 100, nof_repetitions);
comp_times = zeros(1, nof_repetitions);
for i = 1:nof_repetitions
    x = x_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    result = 0; % Pre-allocate memory to store the result (avoids having to find free memory while timing)
    tic
    for j = 1:nof_operations
        % This is the operation being timed
        result = gamma(x);
        % ---------------------------------
    end
    comp_times(i) = toc/nof_operations;
end
comp_times_w_G = comp_times;
fprintf('OK\n')

fprintf('\nEstimating computation time for LEGENDRE FUNCTION ... ')
nof_operations = 1e2;
nof_repetitions = 5e2;
mu_values = round(linspace(-50, 50, nof_repetitions));
v_values = linspace(0, 50, nof_repetitions);
x_values = linspace(0, 100, nof_repetitions);
comp_times = zeros(1, nof_repetitions);
for i = 1:nof_repetitions
    mu = mu_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    v = v_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    x = x_values(i); % Pre-assign value to avoid memory access operation (takes longer)
    result = 0; % Pre-allocate memory to store the result (avoids having to find free memory while timing)
    tic
    for j = 1:nof_operations
        % This is the operation being timed
        result = assoc_legendre_P(mu, v, x, 'double');
        % ---------------------------------
    end
    comp_times(i) = toc/nof_operations;
end
comp_times_w_L = comp_times;
fprintf('OK\n')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear i j x_values y_values mu_values v_values x y mu v result nof_operations nof_repetitions comp_times

% Save results
save('Computation_times.mat')
