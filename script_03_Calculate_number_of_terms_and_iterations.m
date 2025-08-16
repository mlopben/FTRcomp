
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  CALCULATE THE NUMBER OF TERMS AND ITERATIONS FOR VARIOUS TRUNCATION ERRORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% Configuration parameters
m = 5;
K = 5;
Delta = 0.5;
target_trunc_error_values = logspace(-1, -9, 8*4+1);
show_log = 0;

RFM_j_T_values = zeros(size(target_trunc_error_values));
CEM_j_T_values = zeros(size(target_trunc_error_values));
CEM_n_T_values = zeros(size(target_trunc_error_values));
N_iter_RFM_values       = zeros(size(target_trunc_error_values));
N_iter_CEM_1_values     = zeros(size(target_trunc_error_values));
N_iter_CEM_2_num_values = zeros(size(target_trunc_error_values));
N_iter_CEM_2_anl_values = zeros(size(target_trunc_error_values));
N_iter_CEM_3_num_values = zeros(size(target_trunc_error_values));
N_iter_CEM_3_anl_values = zeros(size(target_trunc_error_values));

fprintf('\n                +-----------------------+-----------------------------------------------------------+')
fprintf('\n                |   Config parameters   |                    Number of iteration                    |')
fprintf('\n                +-------+---------------+---------+---------+-------------------+-------------------+')
fprintf('\n                |  RFM  |      CEM      |   RFM   |  CEM-1  |       CEM-2       |       CEM-3       |')
fprintf('\n  +-------------+-------+-------+-------+---------+---------+---------+---------+---------+---------+')
fprintf('\n  |  Trunc Err  |  j_T  |  j_T  |  n_T  |  Ex/An  |  Ex/An  |  Exact  |  Anlyt  |  Exact  |  Anlyt  |')
fprintf('\n  +-------------+-------+-------+-------+---------+---------+---------+---------+---------+---------+')

for i = 1:numel(target_trunc_error_values)
    
    target_trunc_error = target_trunc_error_values(i);
    
    [RFM_j_T, CEM_j_T, CEM_n_T] = calculate_number_of_terms(m, K, Delta, target_trunc_error, show_log);
    
    RFM_j_T_values(i) = RFM_j_T;
    CEM_j_T_values(i) = CEM_j_T;
    CEM_n_T_values(i) = CEM_n_T;
    
    N_iter_RFM       = calculate_number_of_iterations_num(RFM_j_T, -Inf,    m, 'RFM');
    N_iter_CEM_1     = calculate_number_of_iterations_num(CEM_j_T, CEM_n_T, m, 'CEM-I');
    N_iter_CEM_2_num = calculate_number_of_iterations_num(CEM_j_T, CEM_n_T, m, 'CEM-II');
    N_iter_CEM_2_anl = calculate_number_of_iterations_anl(CEM_j_T, CEM_n_T, m, 'CEM-II');
    N_iter_CEM_3_num = calculate_number_of_iterations_num(CEM_j_T, CEM_n_T, m, 'CEM-III');
    N_iter_CEM_3_anl = calculate_number_of_iterations_anl(CEM_j_T, CEM_n_T, m, 'CEM-III');

    N_iter_RFM_values(i)       = N_iter_RFM;
    N_iter_CEM_1_values(i)     = N_iter_CEM_1;
    N_iter_CEM_2_num_values(i) = N_iter_CEM_2_num;
    N_iter_CEM_2_anl_values(i) = N_iter_CEM_2_anl;
    N_iter_CEM_3_num_values(i) = N_iter_CEM_3_num;
    N_iter_CEM_3_anl_values(i) = N_iter_CEM_3_anl;
    
    fprintf('\n  |   %.2e  |   %d  |   %d  |   %d  | %d | %d | %d | %.1f | %d | %.1f |', ...
            target_trunc_error, RFM_j_T, CEM_j_T, CEM_n_T, ...
            N_iter_RFM, N_iter_CEM_1, N_iter_CEM_2_num, N_iter_CEM_2_anl, N_iter_CEM_3_num, N_iter_CEM_3_anl)
    fprintf('\n  +-------------+-------+-------+-------+---------+---------+---------+---------+---------+---------+')

end
fprintf('\n\n')

clear i target_trunc_error RFM_j_T CEM_j_T CEM_n_T N_iter_RFM N_iter_CEM_1 N_iter_CEM_2_num N_iter_CEM_2_anl N_iter_CEM_3_num N_iter_CEM_3_anl

save('Number_of_terms_and_iterations.mat')
