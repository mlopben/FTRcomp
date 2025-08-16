%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IMPACT OF FADING PARAMETER m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1 - Calculate the required number of terms (j_T / n_T)

clear

% Configuration parameters
target_trunc_error = 1e-2;
show_log = 0;

% Default fading parameters
% m = 5;
% K = 5;
% Delta = 0.5;

% Fading parameters
nof_values = 50;
m_values = linspace(1, 20, nof_values);
K = 5;
Delta = 0.5;

RFM_j_T_values = zeros(1, nof_values);
CEM_j_T_values = zeros(1, nof_values);
CEM_n_T_values = zeros(1, nof_values);

for i = 1:nof_values
    m = m_values(i);
    [RFM_j_T_values(i), CEM_j_T_values(i), CEM_n_T_values(i)] = calculate_number_of_terms(m, K, Delta, target_trunc_error, show_log);
    fprintf('m = %.2f\t| RFM_j_T = %d\t| CEM_j_T = %d\t | CEM_n_T = %d\n', m, RFM_j_T_values(i), CEM_j_T_values(i), CEM_n_T_values(i))
end

save('Number_of_terms_for_m_values.mat')

% Step 2 - Calculate computation times based on framework

clear

load('Number_of_terms_for_m_values.mat')

load('Computation_times [Intel Core i9-7900X @ 3.30 GHz][Matlab R2017a].mat', ...
     'comp_times_w_a', 'comp_times_w_m', 'comp_times_w_d', 'comp_times_w_e', 'comp_times_w_r', ...
     'comp_times_w_f', 'comp_times_w_b', 'comp_times_w_G', 'comp_times_w_L')

% Average operation times with median
w_a = median(comp_times_w_a);
w_m = median(comp_times_w_m);
w_d = median(comp_times_w_d);
w_e = median(comp_times_w_e);
w_r = median(comp_times_w_r);
w_f = median(comp_times_w_f);
w_b = median(comp_times_w_b);
w_G = median(comp_times_w_G);
w_L = median(comp_times_w_L);

% Obtain computation times according to the proposed framework
FrameWork_comp_time_RFM_values   = zeros(1, nof_values);
FrameWork_comp_time_CEM_1_values = zeros(1, nof_values);
FrameWork_comp_time_CEM_2_values = zeros(1, nof_values);
FrameWork_comp_time_CEM_3_values = zeros(1, nof_values);

for j = 1:nof_values

    % RFM
    RFM_j_T = RFM_j_T_values(j);
    N_iter_all_j = calculate_number_of_iterations_per_level_anl('N/A', 'j', RFM_j_T, NaN, NaN);
    N_iter_all_k = calculate_number_of_iterations_per_level_anl('N/A', 'k', RFM_j_T, NaN, NaN);
    N_iter_all_l = calculate_number_of_iterations_per_level_anl('N/A', 'l', RFM_j_T, NaN, NaN);
    F_j_all = N_iter_all_j * (w_a);
    F_k_all = N_iter_all_k * (w_a + 2*w_m + w_d + w_e + w_b);
    F_l_RFM = N_iter_all_l * (14*w_a + 8*w_m + 3*w_d + 6*w_e + w_r + w_b + w_G + w_L);
    FrameWork_comp_time_RFM_values(j) = F_j_all + F_k_all + F_l_RFM;
    
    % CEM
    CEM_j_T = CEM_j_T_values(j);
    CEM_n_T = CEM_n_T_values(j);
    N_iter_all_j = calculate_number_of_iterations_per_level_anl('N/A', 'j', CEM_j_T, CEM_n_T, NaN);
    N_iter_all_k = calculate_number_of_iterations_per_level_anl('N/A', 'k', CEM_j_T, CEM_n_T, NaN);
    N_iter_all_l = calculate_number_of_iterations_per_level_anl('N/A', 'l', CEM_j_T, CEM_n_T, NaN);
    N_iter_I_n   = calculate_number_of_iterations_per_level_anl('CEM-I', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_II_n  = calculate_number_of_iterations_per_level_anl('CEM-II', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_III_n = calculate_number_of_iterations_per_level_anl('CEM-III', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_III_p = calculate_number_of_iterations_per_level_anl('CEM-III', 'p', CEM_j_T, CEM_n_T, m);
    F_j_all  = N_iter_all_j * (w_a);
    F_k_all  = N_iter_all_k * (w_a + 2*w_m + w_d + w_e + w_b);
    F_l_CEM1 = N_iter_all_l * (w_a + w_m + w_b);
    F_l_CEM2 = N_iter_all_l * (2*w_a + 2*w_m + w_b);
    F_l_CEM3 = N_iter_all_l * (2*w_a + 2*w_m + w_b);
    F_n_CEM1 = N_iter_I_n   * (15*w_a + 11*w_m + 2*w_d + 2*w_e + w_f + 2*w_G);
    F_n_CEM2 = N_iter_II_n  * (15*w_a + 11*w_m + 2*w_d + 2*w_e + w_f + 2*w_G);
    F_n_CEM3 = N_iter_III_n * (11*w_a + 7*w_m + 2*w_d + 2*w_e + w_f);
    F_p_CEM3 = N_iter_III_p * (3*w_a + 2*w_m);
    FrameWork_comp_time_CEM_1_values(j) = F_j_all + F_k_all + F_l_CEM1 + F_n_CEM1;
    FrameWork_comp_time_CEM_2_values(j) = F_j_all + F_k_all + F_l_CEM2 + F_n_CEM2;
    FrameWork_comp_time_CEM_3_values(j) = F_j_all + F_k_all + F_l_CEM3 + F_n_CEM3 + F_p_CEM3;
    
end

% Step 3 - Plot results
FontSize = 14;
figure
plot(m_values, FrameWork_comp_time_RFM_values, 'k-', ...
     m_values, FrameWork_comp_time_CEM_1_values, 'r-', ...
     m_values, FrameWork_comp_time_CEM_2_values, 'g-', ...
     m_values, FrameWork_comp_time_CEM_3_values, 'b-')
set(gca, 'YScale', 'log')
box on
grid
set(gca, 'YTick', 10.^(-2:2))
set(gca, 'YLim', [10^-2 10^2])
set(gca, 'XTick', 0:2:20)
set(gca, 'XLim', [0 20])
set(gca, 'FontSize', FontSize)
xlabel('Fading parameter $m$', 'Interpreter', 'LaTeX')
ylabel('Computation time (s)', 'Interpreter', 'LaTeX')
title('Intel Core i9-7900X @ 3.30 GHz', 'Interpreter', 'LaTeX')
h = legend('RFM', 'CEM-I', 'CEM-II', 'CEM-III');
set(h, 'Location', 'East', 'Interpreter', 'LaTeX', 'FontSize', FontSize-4);

% Save the figure
print('-depsc', '-tiff', '-r600', 'Computation_times_vs_m.eps')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print('-dpdf', 'Computation_times_vs_m.pdf')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IMPACT OF FADING PARAMETER K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1 - Calculate the required number of terms (j_T / n_T)

clear

% Configuration parameters
target_trunc_error = 1e-2;
show_log = 0;

% Default fading parameters
% m = 5;
% K = 5;
% Delta = 0.5;

% Fading parameters
nof_values = 50;
m = 5;
K_values = linspace(0.02, 10, nof_values);
Delta = 0.5;

RFM_j_T_values = zeros(1, nof_values);
CEM_j_T_values = zeros(1, nof_values);
CEM_n_T_values = zeros(1, nof_values);

for i = 1:nof_values
    K = K_values(i);
    [RFM_j_T_values(i), CEM_j_T_values(i), CEM_n_T_values(i)] = calculate_number_of_terms(m, K, Delta, target_trunc_error, show_log);
    fprintf('K = %.2f\t| RFM_j_T = %d\t| CEM_j_T = %d\t | CEM_n_T = %d\n', K, RFM_j_T_values(i), CEM_j_T_values(i), CEM_n_T_values(i))
end

save('Number_of_terms_for_K_values.mat')

% Step 2 - Calculate computation times based on framework

clear

load('Number_of_terms_for_K_values.mat')

load('Computation_times [Intel Core i9-7900X @ 3.30 GHz][Matlab R2017a].mat', ...
     'comp_times_w_a', 'comp_times_w_m', 'comp_times_w_d', 'comp_times_w_e', 'comp_times_w_r', ...
     'comp_times_w_f', 'comp_times_w_b', 'comp_times_w_G', 'comp_times_w_L')

% Average operation times with median
w_a = median(comp_times_w_a);
w_m = median(comp_times_w_m);
w_d = median(comp_times_w_d);
w_e = median(comp_times_w_e);
w_r = median(comp_times_w_r);
w_f = median(comp_times_w_f);
w_b = median(comp_times_w_b);
w_G = median(comp_times_w_G);
w_L = median(comp_times_w_L);

% Obtain computation times according to the proposed framework
FrameWork_comp_time_RFM_values   = zeros(1, nof_values);
FrameWork_comp_time_CEM_1_values = zeros(1, nof_values);
FrameWork_comp_time_CEM_2_values = zeros(1, nof_values);
FrameWork_comp_time_CEM_3_values = zeros(1, nof_values);

for j = 1:nof_values

    % RFM
    RFM_j_T = RFM_j_T_values(j);
    N_iter_all_j = calculate_number_of_iterations_per_level_anl('N/A', 'j', RFM_j_T, NaN, NaN);
    N_iter_all_k = calculate_number_of_iterations_per_level_anl('N/A', 'k', RFM_j_T, NaN, NaN);
    N_iter_all_l = calculate_number_of_iterations_per_level_anl('N/A', 'l', RFM_j_T, NaN, NaN);
    F_j_all = N_iter_all_j * (w_a);
    F_k_all = N_iter_all_k * (w_a + 2*w_m + w_d + w_e + w_b);
    F_l_RFM = N_iter_all_l * (14*w_a + 8*w_m + 3*w_d + 6*w_e + w_r + w_b + w_G + w_L);
    FrameWork_comp_time_RFM_values(j) = F_j_all + F_k_all + F_l_RFM;
    
    % CEM
    CEM_j_T = CEM_j_T_values(j);
    CEM_n_T = CEM_n_T_values(j);
    N_iter_all_j = calculate_number_of_iterations_per_level_anl('N/A', 'j', CEM_j_T, CEM_n_T, NaN);
    N_iter_all_k = calculate_number_of_iterations_per_level_anl('N/A', 'k', CEM_j_T, CEM_n_T, NaN);
    N_iter_all_l = calculate_number_of_iterations_per_level_anl('N/A', 'l', CEM_j_T, CEM_n_T, NaN);
    N_iter_I_n   = calculate_number_of_iterations_per_level_anl('CEM-I', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_II_n  = calculate_number_of_iterations_per_level_anl('CEM-II', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_III_n = calculate_number_of_iterations_per_level_anl('CEM-III', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_III_p = calculate_number_of_iterations_per_level_anl('CEM-III', 'p', CEM_j_T, CEM_n_T, m);
    F_j_all  = N_iter_all_j * (w_a);
    F_k_all  = N_iter_all_k * (w_a + 2*w_m + w_d + w_e + w_b);
    F_l_CEM1 = N_iter_all_l * (w_a + w_m + w_b);
    F_l_CEM2 = N_iter_all_l * (2*w_a + 2*w_m + w_b);
    F_l_CEM3 = N_iter_all_l * (2*w_a + 2*w_m + w_b);
    F_n_CEM1 = N_iter_I_n   * (15*w_a + 11*w_m + 2*w_d + 2*w_e + w_f + 2*w_G);
    F_n_CEM2 = N_iter_II_n  * (15*w_a + 11*w_m + 2*w_d + 2*w_e + w_f + 2*w_G);
    F_n_CEM3 = N_iter_III_n * (11*w_a + 7*w_m + 2*w_d + 2*w_e + w_f);
    F_p_CEM3 = N_iter_III_p * (3*w_a + 2*w_m);
    FrameWork_comp_time_CEM_1_values(j) = F_j_all + F_k_all + F_l_CEM1 + F_n_CEM1;
    FrameWork_comp_time_CEM_2_values(j) = F_j_all + F_k_all + F_l_CEM2 + F_n_CEM2;
    FrameWork_comp_time_CEM_3_values(j) = F_j_all + F_k_all + F_l_CEM3 + F_n_CEM3 + F_p_CEM3;
    
end

% Step 3 - Plot results
FontSize = 14;
figure
plot(K_values, FrameWork_comp_time_RFM_values, 'k-', ...
     K_values, FrameWork_comp_time_CEM_1_values, 'r-', ...
     K_values, FrameWork_comp_time_CEM_2_values, 'g-', ...
     K_values, FrameWork_comp_time_CEM_3_values, 'b-')
set(gca, 'YScale', 'log')
box on
grid
set(gca, 'YTick', 10.^(-4:2))
set(gca, 'YLim', [10^-4 10^2])
set(gca, 'XTick', 0:1:10)
set(gca, 'XLim', [0 10])
set(gca, 'FontSize', FontSize)
xlabel('Fading parameter $K$', 'Interpreter', 'LaTeX')
ylabel('Computation time (s)', 'Interpreter', 'LaTeX')
title('Intel Core i9-7900X @ 3.30 GHz', 'Interpreter', 'LaTeX')
h = legend('RFM', 'CEM-I', 'CEM-II', 'CEM-III');
set(h, 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'FontSize', FontSize-4);

% Save the figure
print('-depsc', '-tiff', '-r600', 'Computation_times_vs_K.eps')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print('-dpdf', 'Computation_times_vs_K.pdf')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  IMPACT OF FADING PARAMETER Delta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1 - Calculate the required number of terms (j_T / n_T)

clear

% Configuration parameters
target_trunc_error = 1e-2;
show_log = 0;

% Default fading parameters
% m = 5;
% K = 5;
% Delta = 0.5;

% Fading parameters
nof_values = 50;
m = 5;
K = 5;
Delta_values = linspace(0.001, 1, nof_values);

RFM_j_T_values = zeros(1, nof_values);
CEM_j_T_values = zeros(1, nof_values);
CEM_n_T_values = zeros(1, nof_values);

for i = 1:nof_values
    Delta = Delta_values(i);
    [RFM_j_T_values(i), CEM_j_T_values(i), CEM_n_T_values(i)] = calculate_number_of_terms(m, K, Delta, target_trunc_error, show_log);
    fprintf('Delta = %.2f\t| RFM_j_T = %d\t| CEM_j_T = %d\t | CEM_n_T = %d\n', Delta, RFM_j_T_values(i), CEM_j_T_values(i), CEM_n_T_values(i))
end

save('Number_of_terms_for_Delta_values.mat')

% Step 2 - Calculate computation times based on framework

clear

load('Number_of_terms_for_Delta_values.mat')

load('Computation_times [Intel Core i9-7900X @ 3.30 GHz][Matlab R2017a].mat', ...
     'comp_times_w_a', 'comp_times_w_m', 'comp_times_w_d', 'comp_times_w_e', 'comp_times_w_r', ...
     'comp_times_w_f', 'comp_times_w_b', 'comp_times_w_G', 'comp_times_w_L')

% Average operation times with median
w_a = median(comp_times_w_a);
w_m = median(comp_times_w_m);
w_d = median(comp_times_w_d);
w_e = median(comp_times_w_e);
w_r = median(comp_times_w_r);
w_f = median(comp_times_w_f);
w_b = median(comp_times_w_b);
w_G = median(comp_times_w_G);
w_L = median(comp_times_w_L);

% Obtain computation times according to the proposed framework
FrameWork_comp_time_RFM_values   = zeros(1, nof_values);
FrameWork_comp_time_CEM_1_values = zeros(1, nof_values);
FrameWork_comp_time_CEM_2_values = zeros(1, nof_values);
FrameWork_comp_time_CEM_3_values = zeros(1, nof_values);

for j = 1:nof_values

    % RFM
    RFM_j_T = RFM_j_T_values(j);
    N_iter_all_j = calculate_number_of_iterations_per_level_anl('N/A', 'j', RFM_j_T, NaN, NaN);
    N_iter_all_k = calculate_number_of_iterations_per_level_anl('N/A', 'k', RFM_j_T, NaN, NaN);
    N_iter_all_l = calculate_number_of_iterations_per_level_anl('N/A', 'l', RFM_j_T, NaN, NaN);
    F_j_all = N_iter_all_j * (w_a);
    F_k_all = N_iter_all_k * (w_a + 2*w_m + w_d + w_e + w_b);
    F_l_RFM = N_iter_all_l * (14*w_a + 8*w_m + 3*w_d + 6*w_e + w_r + w_b + w_G + w_L);
    FrameWork_comp_time_RFM_values(j) = F_j_all + F_k_all + F_l_RFM;
    
    % CEM
    CEM_j_T = CEM_j_T_values(j);
    CEM_n_T = CEM_n_T_values(j);
    N_iter_all_j = calculate_number_of_iterations_per_level_anl('N/A', 'j', CEM_j_T, CEM_n_T, NaN);
    N_iter_all_k = calculate_number_of_iterations_per_level_anl('N/A', 'k', CEM_j_T, CEM_n_T, NaN);
    N_iter_all_l = calculate_number_of_iterations_per_level_anl('N/A', 'l', CEM_j_T, CEM_n_T, NaN);
    N_iter_I_n   = calculate_number_of_iterations_per_level_anl('CEM-I', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_II_n  = calculate_number_of_iterations_per_level_anl('CEM-II', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_III_n = calculate_number_of_iterations_per_level_anl('CEM-III', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_III_p = calculate_number_of_iterations_per_level_anl('CEM-III', 'p', CEM_j_T, CEM_n_T, m);
    F_j_all  = N_iter_all_j * (w_a);
    F_k_all  = N_iter_all_k * (w_a + 2*w_m + w_d + w_e + w_b);
    F_l_CEM1 = N_iter_all_l * (w_a + w_m + w_b);
    F_l_CEM2 = N_iter_all_l * (2*w_a + 2*w_m + w_b);
    F_l_CEM3 = N_iter_all_l * (2*w_a + 2*w_m + w_b);
    F_n_CEM1 = N_iter_I_n   * (15*w_a + 11*w_m + 2*w_d + 2*w_e + w_f + 2*w_G);
    F_n_CEM2 = N_iter_II_n  * (15*w_a + 11*w_m + 2*w_d + 2*w_e + w_f + 2*w_G);
    F_n_CEM3 = N_iter_III_n * (11*w_a + 7*w_m + 2*w_d + 2*w_e + w_f);
    F_p_CEM3 = N_iter_III_p * (3*w_a + 2*w_m);
    FrameWork_comp_time_CEM_1_values(j) = F_j_all + F_k_all + F_l_CEM1 + F_n_CEM1;
    FrameWork_comp_time_CEM_2_values(j) = F_j_all + F_k_all + F_l_CEM2 + F_n_CEM2;
    FrameWork_comp_time_CEM_3_values(j) = F_j_all + F_k_all + F_l_CEM3 + F_n_CEM3 + F_p_CEM3;
    
end

% Step 3 - Plot results
FontSize = 14;
figure
plot(Delta_values, FrameWork_comp_time_RFM_values, 'k-', ...
     Delta_values, FrameWork_comp_time_CEM_1_values, 'r-', ...
     Delta_values, FrameWork_comp_time_CEM_2_values, 'g-', ...
     Delta_values, FrameWork_comp_time_CEM_3_values, 'b-')
set(gca, 'YScale', 'log')
box on
grid
set(gca, 'YTick', 10.^(-3:2))
set(gca, 'YLim', [10^-3 10^2])
set(gca, 'XTick', 0:0.1:1)
set(gca, 'XLim', [0 1])
set(gca, 'FontSize', FontSize)
xlabel('Fading parameter $\Delta$', 'Interpreter', 'LaTeX')
ylabel('Computation time (s)', 'Interpreter', 'LaTeX')
title('Intel Core i9-7900X @ 3.30 GHz', 'Interpreter', 'LaTeX')
h = legend('RFM', 'CEM-I', 'CEM-II', 'CEM-III');
set(h, 'Location', 'SouthEast', 'Interpreter', 'LaTeX', 'FontSize', FontSize-4);

% Save the figure
print('-depsc', '-tiff', '-r600', 'Computation_times_vs_Delta.eps')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print('-dpdf', 'Computation_times_vs_Delta.pdf')
