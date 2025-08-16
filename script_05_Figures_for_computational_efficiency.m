
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of iterations vs. truncation error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

FontSize = 14;

load('Number_of_terms_and_iterations.mat')

fprintf('\nMax relative error of N_iter_II_n, N_iter_III_n [equation (17)] = %.2e\n', max(abs(N_iter_CEM_2_num_values - N_iter_CEM_2_anl_values)./N_iter_CEM_2_num_values))
fprintf('\nMax relative error of N_iter_III_p [equation (24)] = %.2e\n', max(abs(N_iter_CEM_3_num_values - N_iter_CEM_3_anl_values)./N_iter_CEM_3_num_values))

N_iter_R_j   = zeros(size(target_trunc_error_values));
N_iter_R_k   = zeros(size(target_trunc_error_values));
N_iter_R_l   = zeros(size(target_trunc_error_values));
N_iter_C_j   = zeros(size(target_trunc_error_values));
N_iter_C_k   = zeros(size(target_trunc_error_values));
N_iter_C_l   = zeros(size(target_trunc_error_values));
N_iter_I_n   = zeros(size(target_trunc_error_values));
N_iter_II_n  = zeros(size(target_trunc_error_values));
N_iter_III_n = zeros(size(target_trunc_error_values));
N_iter_III_p = zeros(size(target_trunc_error_values));

for i = 1:numel(target_trunc_error_values)

    % RFM
    RFM_j_T = RFM_j_T_values(i);
    N_iter_R_j(i) = calculate_number_of_iterations_per_level_anl('N/A', 'j', RFM_j_T, NaN, NaN);
    N_iter_R_k(i) = calculate_number_of_iterations_per_level_anl('N/A', 'k', RFM_j_T, NaN, NaN);
    N_iter_R_l(i) = calculate_number_of_iterations_per_level_anl('N/A', 'l', RFM_j_T, NaN, NaN);

    % CEM
    CEM_j_T = CEM_j_T_values(i);
    CEM_n_T = CEM_n_T_values(i);
    N_iter_C_j(i)   = calculate_number_of_iterations_per_level_anl('N/A', 'j', CEM_j_T, CEM_n_T, NaN);
    N_iter_C_k(i)   = calculate_number_of_iterations_per_level_anl('N/A', 'k', CEM_j_T, CEM_n_T, NaN);
    N_iter_C_l(i)   = calculate_number_of_iterations_per_level_anl('N/A', 'l', CEM_j_T, CEM_n_T, NaN);
    N_iter_I_n(i)   = calculate_number_of_iterations_per_level_anl('CEM-I', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_II_n(i)  = calculate_number_of_iterations_per_level_anl('CEM-II', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_III_n(i) = calculate_number_of_iterations_per_level_anl('CEM-III', 'n', CEM_j_T, CEM_n_T, NaN);
    N_iter_III_p(i) = calculate_number_of_iterations_per_level_anl('CEM-III', 'p', CEM_j_T, CEM_n_T, m);

end

% Notice that: N_iter_R_l   == N_iter_RFM_values
% Notice that: N_iter_I_n   == N_iter_CEM_1_values
% Notice that: N_iter_II_n  == N_iter_CEM_2_anl_values
% Notice that: N_iter_III_p == N_iter_CEM_3_anl_values
figure
h = loglog(target_trunc_error_values, N_iter_R_j, 'k-', ...
           target_trunc_error_values, N_iter_R_k, 'k--', ...
           target_trunc_error_values, N_iter_R_l, 'k-.', ...                % == N_iter_RFM_values
           target_trunc_error_values, N_iter_C_j, 'r-', ...
           target_trunc_error_values, N_iter_C_k, 'r--', ...
           target_trunc_error_values, N_iter_C_l, 'r-.', ...
           target_trunc_error_values, N_iter_I_n, 'b-', ...                 % == N_iter_CEM_1_values
           target_trunc_error_values, N_iter_CEM_2_num_values, 'b--', ...
           target_trunc_error_values, N_iter_II_n, 'bo', ...                % == N_iter_CEM_2_anl_values
           target_trunc_error_values, N_iter_CEM_3_num_values, 'm-', ...
           target_trunc_error_values, N_iter_III_p, 'ms');                  % == N_iter_CEM_3_anl_values
set(h(9), 'MarkerSize', 4)
set(h(11), 'MarkerSize', 5)
box on
grid
set(gca, 'YTick', 10.^(0:8))
set(gca, 'YLim', [10^0 10^8])
set(gca, 'XTick', 10.^(-9:-1))
set(gca, 'XLim', [target_trunc_error_values(end) target_trunc_error_values(1)])
set(gca, 'FontSize', FontSize)
xlabel('Maximum truncation error, $\varepsilon_T^{\max}$', 'Interpreter', 'LaTeX')
ylabel('Number of iterations, $N_\mathrm{iter}$', 'Interpreter', 'LaTeX')

h = legend('$N_\mathrm{iter}^{\mathrm{all},j}$ (RFM)', ...
           '$N_\mathrm{iter}^{\mathrm{all},k}$ (RFM)', ...
           '$N_\mathrm{iter}^{\mathrm{all},l}$ (RFM)', ...
           '$N_\mathrm{iter}^{\mathrm{all},j}$ (CEMs)', ...
           '$N_\mathrm{iter}^{\mathrm{all},k}$ (CEMs)', ...
           '$N_\mathrm{iter}^{\mathrm{all},l}$ (CEMs)', ...
           '$N_\mathrm{iter}^{\mathrm{I},n}$', ...
           '$N_\mathrm{iter}^{\mathrm{II},n}, N_\mathrm{iter}^{\mathrm{III},n}$ (Exact)', ...
           '$N_\mathrm{iter}^{\mathrm{II},n}, N_\mathrm{iter}^{\mathrm{III},n}$ (Anl.)', ...
           '$N_\mathrm{iter}^{\mathrm{III},p}$ (Exact)', ...
           '$N_\mathrm{iter}^{\mathrm{III},p}$ (Anl.)');
set(h, 'Location', 'SouthWest', 'Interpreter', 'LaTeX', 'FontSize', FontSize-4);

% Save the figure
print -depsc -tiff -r600 Number_of_iterations.eps
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print -dpdf Number_of_iterations.pdf


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation times vs. truncation error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

% FontSize = 18;
FontSize = 14;

processors = {'Intel Core i3-3220 @ 3.30 GHz';
              'Intel Core i5-6500 @ 3.20 GHz';
              'Intel Core i7-2600 @ 3.40 GHz';
              'Intel Core i7-8700 @ 3.20 GHz';
              'Intel Core i9-7900X @ 3.30 GHz'};

data_files = {'Computation_times [Intel Core i3-3220 @ 3.30 GHz][Matlab R2017a].mat';
              'Computation_times [Intel Core i5-6500 @ 3.20 GHz][Matlab R2017a].mat';
              'Computation_times [Intel Core i7-2600 @ 3.40 GHz][Matlab R2017a].mat';
              'Computation_times [Intel Core i7-8700 @ 3.20 GHz][Matlab R2017a].mat';
              'Computation_times [Intel Core i9-7900X @ 3.30 GHz][Matlab R2017a].mat'};

% Select how the final value is selected based on the different repetitions
% experimental_time_calculation = 'mean';
experimental_time_calculation = 'median';
% experimental_time_calculation = 'min';

for i = 1:numel(data_files)
    
    % Load data file
    load(data_files{i})
    
    % Plot experimental computation times for each fading model
    figure(i)
    switch experimental_time_calculation
        case 'mean'
            error('Use median calcualtion')
            loglog(target_trunc_error_values, mean(comp_time_RFM_values,1), 'k-', ...
                   target_trunc_error_values, mean(comp_time_CEM_1_values,1), 'r-', ...
                   target_trunc_error_values, mean(comp_time_CEM_2_values,1), 'g-', ...
                   target_trunc_error_values, mean(comp_time_CEM_3_values,1), 'b-')
        case 'median'
            loglog(target_trunc_error_values, median(comp_time_RFM_values,1), 'k-', ...
                   target_trunc_error_values, median(comp_time_CEM_1_values,1), 'r-', ...
                   target_trunc_error_values, median(comp_time_CEM_2_values,1), 'g-', ...
                   target_trunc_error_values, median(comp_time_CEM_3_values,1), 'b-')
        case 'min'
            error('Use median calcualtion')
            loglog(target_trunc_error_values, min(comp_time_RFM_values,[],1), 'k-', ...
                   target_trunc_error_values, min(comp_time_CEM_1_values,[],1), 'r-', ...
                   target_trunc_error_values, min(comp_time_CEM_2_values,[],1), 'g-', ...
                   target_trunc_error_values, min(comp_time_CEM_3_values,[],1), 'b-')
        otherwise
            error('Wrong calculation type')
    end
    box on
    grid
    set(gca, 'YTick', 10.^(-3:3))
    set(gca, 'YLim', [10^-3 10^3])
    set(gca, 'XTick', 10.^(-9:-1))
    set(gca, 'XLim', [target_trunc_error_values(end) target_trunc_error_values(1)])
    set(gca, 'FontSize', FontSize)
    xlabel('Maximum truncation error, $\varepsilon_T^{\max}$', 'Interpreter', 'LaTeX')
    ylabel('Computation time (s)', 'Interpreter', 'LaTeX')
    title(processors{i}, 'Interpreter', 'LaTeX')
    
    % % Plot FLOP times for each operation
    % figure;
    % subplot(3,3,1); plot(comp_times_w_a); title(sprintf('w_a: %.2es | %.2es | %.2es', mean(comp_times_w_a), median(comp_times_w_a), min(comp_times_w_a)))
    % subplot(3,3,2); plot(comp_times_w_m); title(sprintf('w_m: %.2es | %.2es | %.2es', mean(comp_times_w_m), median(comp_times_w_m), min(comp_times_w_m)))
    % subplot(3,3,3); plot(comp_times_w_d); title(sprintf('w_d: %.2es | %.2es | %.2es', mean(comp_times_w_d), median(comp_times_w_d), min(comp_times_w_d)))
    % subplot(3,3,4); plot(comp_times_w_e); title(sprintf('w_e: %.2es | %.2es | %.2es', mean(comp_times_w_e), median(comp_times_w_e), min(comp_times_w_e)))
    % subplot(3,3,5); plot(comp_times_w_r); title(sprintf('w_r: %.2es | %.2es | %.2es', mean(comp_times_w_r), median(comp_times_w_r), min(comp_times_w_r)))
    % subplot(3,3,6); plot(comp_times_w_f); title(sprintf('w_f: %.2es | %.2es | %.2es', mean(comp_times_w_f), median(comp_times_w_f), min(comp_times_w_f)))
    % subplot(3,3,7); plot(comp_times_w_b); title(sprintf('w_b: %.2es | %.2es | %.2es', mean(comp_times_w_b), median(comp_times_w_b), min(comp_times_w_b)))
    % subplot(3,3,8); plot(comp_times_w_G); title(sprintf('w_G: %.2es | %.2es | %.2es', mean(comp_times_w_G), median(comp_times_w_G), min(comp_times_w_G)))
    % subplot(3,3,9); plot(comp_times_w_L); title(sprintf('w_L: %.2es | %.2es | %.2es', mean(comp_times_w_L), median(comp_times_w_L), min(comp_times_w_L)))
    % fprintf('\n\nFLOP times - %s', processors{i})
    % fprintf('\nSave the figure then press a key to continue...\n\n')
    % pause
    % close(gcf)
    
    % Get average (mean/median) FLOP times for each operation
    fprintf('\nFLOP times (seconds) for %s:', processors{i})
    fprintf('\n------------------------------------------------------------------------')
    fprintf('\nw_a: %.2e [%.2e/%d] \t\t\t| %.2e [%.2e/%d] \t\t\t| %.2e [%.2e/%d]', mean(comp_times_w_a), mean(comp_times_w_a)/mean(comp_times_w_a), round(mean(comp_times_w_a)/mean(comp_times_w_a)), median(comp_times_w_a), median(comp_times_w_a)/median(comp_times_w_a), round(median(comp_times_w_a)/median(comp_times_w_a)), min(comp_times_w_a), min(comp_times_w_a)/min(comp_times_w_a), round(min(comp_times_w_a)/min(comp_times_w_a)))
    fprintf('\nw_m: %.2e [%.2e/%d] \t\t\t| %.2e [%.2e/%d] \t\t\t| %.2e [%.2e/%d]', mean(comp_times_w_m), mean(comp_times_w_m)/mean(comp_times_w_a), round(mean(comp_times_w_m)/mean(comp_times_w_a)), median(comp_times_w_m), median(comp_times_w_m)/median(comp_times_w_a), round(median(comp_times_w_m)/median(comp_times_w_a)), min(comp_times_w_m), min(comp_times_w_m)/min(comp_times_w_a), round(min(comp_times_w_m)/min(comp_times_w_a)))
    fprintf('\nw_d: %.2e [%.2e/%d] \t\t\t| %.2e [%.2e/%d] \t\t\t| %.2e [%.2e/%d]', mean(comp_times_w_d), mean(comp_times_w_d)/mean(comp_times_w_a), round(mean(comp_times_w_d)/mean(comp_times_w_a)), median(comp_times_w_d), median(comp_times_w_d)/median(comp_times_w_a), round(median(comp_times_w_d)/median(comp_times_w_a)), min(comp_times_w_d), min(comp_times_w_d)/min(comp_times_w_a), round(min(comp_times_w_d)/min(comp_times_w_a)))
    fprintf('\nw_e: %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d]',   mean(comp_times_w_e), mean(comp_times_w_e)/mean(comp_times_w_a), round(mean(comp_times_w_e)/mean(comp_times_w_a)), median(comp_times_w_e), median(comp_times_w_e)/median(comp_times_w_a), round(median(comp_times_w_e)/median(comp_times_w_a)), min(comp_times_w_e), min(comp_times_w_e)/min(comp_times_w_a), round(min(comp_times_w_e)/min(comp_times_w_a)))
    fprintf('\nw_r: %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d]',   mean(comp_times_w_r), mean(comp_times_w_r)/mean(comp_times_w_a), round(mean(comp_times_w_r)/mean(comp_times_w_a)), median(comp_times_w_r), median(comp_times_w_r)/median(comp_times_w_a), round(median(comp_times_w_r)/median(comp_times_w_a)), min(comp_times_w_r), min(comp_times_w_r)/min(comp_times_w_a), round(min(comp_times_w_r)/min(comp_times_w_a)))
    fprintf('\nw_f: %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d]',   mean(comp_times_w_f), mean(comp_times_w_f)/mean(comp_times_w_a), round(mean(comp_times_w_f)/mean(comp_times_w_a)), median(comp_times_w_f), median(comp_times_w_f)/median(comp_times_w_a), round(median(comp_times_w_f)/median(comp_times_w_a)), min(comp_times_w_f), min(comp_times_w_f)/min(comp_times_w_a), round(min(comp_times_w_f)/min(comp_times_w_a)))
    fprintf('\nw_b: %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d]',   mean(comp_times_w_b), mean(comp_times_w_b)/mean(comp_times_w_a), round(mean(comp_times_w_b)/mean(comp_times_w_a)), median(comp_times_w_b), median(comp_times_w_b)/median(comp_times_w_a), round(median(comp_times_w_b)/median(comp_times_w_a)), min(comp_times_w_b), min(comp_times_w_b)/min(comp_times_w_a), round(min(comp_times_w_b)/min(comp_times_w_a)))
    fprintf('\nw_G: %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d] \t\t| %.2e [%.2e/%d]',   mean(comp_times_w_G), mean(comp_times_w_G)/mean(comp_times_w_a), round(mean(comp_times_w_G)/mean(comp_times_w_a)), median(comp_times_w_G), median(comp_times_w_G)/median(comp_times_w_a), round(median(comp_times_w_G)/median(comp_times_w_a)), min(comp_times_w_G), min(comp_times_w_G)/min(comp_times_w_a), round(min(comp_times_w_G)/min(comp_times_w_a)))
    fprintf('\nw_L: %.2e [%.2e/%d] \t| %.2e [%.2e/%d] \t| %.2e [%.2e/%d]',     mean(comp_times_w_L), mean(comp_times_w_L)/mean(comp_times_w_a), round(mean(comp_times_w_L)/mean(comp_times_w_a)), median(comp_times_w_L), median(comp_times_w_L)/median(comp_times_w_a), round(median(comp_times_w_L)/median(comp_times_w_a)), min(comp_times_w_L), min(comp_times_w_L)/min(comp_times_w_a), round(min(comp_times_w_L)/min(comp_times_w_a)))
    fprintf('\n\n\n')
    switch experimental_time_calculation
        case 'mean'
            error('Use median calcualtion')
            w_a = mean(comp_times_w_a);
            w_m = mean(comp_times_w_m);
            w_d = mean(comp_times_w_d);
            w_e = mean(comp_times_w_e);
            w_r = mean(comp_times_w_r);
            w_f = mean(comp_times_w_f);
            w_b = mean(comp_times_w_b);
            w_G = mean(comp_times_w_G);
            w_L = mean(comp_times_w_L);
        case 'median'
            w_a = median(comp_times_w_a);
            w_m = median(comp_times_w_m);
            w_d = median(comp_times_w_d);
            w_e = median(comp_times_w_e);
            w_r = median(comp_times_w_r);
            w_f = median(comp_times_w_f);
            w_b = median(comp_times_w_b);
            w_G = median(comp_times_w_G);
            w_L = median(comp_times_w_L);
        case 'min'
            error('Use median calcualtion')
            w_a = min(comp_times_w_a);
            w_m = min(comp_times_w_m);
            w_d = min(comp_times_w_d);
            w_e = min(comp_times_w_e);
            w_r = min(comp_times_w_r);
            w_f = min(comp_times_w_f);
            w_b = min(comp_times_w_b);
            w_G = min(comp_times_w_G);
            w_L = min(comp_times_w_L);
        otherwise
            error('Wrong calculation type')
    end
    
    % Obtain computation times according to the proposed framework
    FrameWork_comp_time_RFM_values   = zeros(size(target_trunc_error_values));
    FrameWork_comp_time_CEM_1_values = zeros(size(target_trunc_error_values));
    FrameWork_comp_time_CEM_2_values = zeros(size(target_trunc_error_values));
    FrameWork_comp_time_CEM_3_values = zeros(size(target_trunc_error_values));
    
    for j = 1:numel(target_trunc_error_values)
        
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
    
    % Add results to the figure
    figure(i)
    hold
    plot(target_trunc_error_values, FrameWork_comp_time_RFM_values, 'k--', ...
         target_trunc_error_values, FrameWork_comp_time_CEM_1_values, 'r--', ...
         target_trunc_error_values, FrameWork_comp_time_CEM_2_values, 'g--', ...
         target_trunc_error_values, FrameWork_comp_time_CEM_3_values, 'b--')
    h = legend('RFM (experimental)', 'CEM-I (experimental)', 'CEM-II (experimental)', 'CEM-III (experimental)', ...
               'RFM (framework)', 'CEM-I (framework)', 'CEM-II (framework)', 'CEM-III (framework)');
    set(h, 'Location', 'SouthWest', 'Interpreter', 'LaTeX', 'FontSize', FontSize-4);
    
    % Save the figure
    file_name = processors{i};
    file_name(strfind(file_name, ' ')) = '_';
    file_name(strfind(file_name, '-')) = '_';
    file_name(strfind(file_name, '.')) = '_';
    
    print('-depsc', '-tiff', '-r600', ['Computation_times_' file_name '.eps'])
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
    print('-dpdf', ['Computation_times_' file_name '.pdf'])

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalised FLOP times vs. input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

FontSize = 12;

load('Computation_times [Intel Core i5-6500 @ 3.20 GHz][Matlab R2017a].mat');

subplot(2,2,1)
x_values = round(linspace(0, 15, 500));
h = plot(x_values, ones(size(x_values)), 'k--', ...
         x_values, comp_times_w_f./median(comp_times_w_f), 'b-');
set(h(1), 'Color', [0.3 0.3 0.3])
box on
grid minor
set(gca, 'XLim', [min(x_values) max(x_values)])
set(gca, 'XTick', min(x_values):5:max(x_values))
set(gca, 'YLim', [0.4 1.1])
set(gca, 'YTick', [0.5 1])
set(gca, 'FontSize', FontSize)
xlabel('$x$', 'Interpreter', 'LaTeX')
ylabel('Normalised time', 'Interpreter', 'LaTeX')
title('Factorial, $x!$', 'Interpreter', 'LaTeX')

subplot(2,2,2)
x_values = round(linspace(0, 50, 500));
h = plot(x_values, ones(size(x_values)), 'k--', ...
         x_values, comp_times_w_b./median(comp_times_w_b), 'b-');
set(h(1), 'Color', [0.3 0.3 0.3])
box on
grid minor
set(gca, 'XLim', [min(x_values) max(x_values)])
set(gca, 'XTick', min(x_values):10:max(x_values))
set(gca, 'YLim', [0.2 1.2])
set(gca, 'YTick', [0.5 1])
set(gca, 'FontSize', FontSize)
xlabel('$x$', 'Interpreter', 'LaTeX')
ylabel('Normalised time', 'Interpreter', 'LaTeX')
title('Binomial coefficient, ${x \choose \, x/2 \, }$', 'Interpreter', 'LaTeX')

subplot(2,2,3)
x_values = linspace(0, 100, 500);
h = plot(x_values, ones(size(x_values)), 'k--', ...
         x_values, comp_times_w_G./median(comp_times_w_G), 'b-');
set(h(1), 'Color', [0.3 0.3 0.3])
box on
grid minor
set(gca, 'XLim', [min(x_values) max(x_values)])
set(gca, 'XTick', min(x_values):20:max(x_values))
set(gca, 'YLim', [0.8 1.1])
set(gca, 'YTick', [0.9 1])
set(gca, 'FontSize', FontSize)
xlabel('$x$', 'Interpreter', 'LaTeX')
ylabel('Normalised time', 'Interpreter', 'LaTeX')
title('Gamma function, $\Gamma(x)$', 'Interpreter', 'LaTeX')

subplot(2,2,4)
mu_values = round(linspace(-50, 50, 500));
h = plot(mu_values, ones(size(mu_values)), 'k--', ...
         mu_values, comp_times_w_L./median(comp_times_w_L), 'b-');
set(h(1), 'Color', [0.3 0.3 0.3])
box on
grid minor
set(gca, 'XLim', [min(mu_values) max(mu_values)])
set(gca, 'XTick', min(mu_values):20:max(mu_values))
set(gca, 'YLim', [0 3])
set(gca, 'YTick', 0:3)
set(gca, 'FontSize', FontSize)
xlabel('$\mu$', 'Interpreter', 'LaTeX')
ylabel('Normalised time', 'Interpreter', 'LaTeX')
title('Legendre function, $P_\nu^\mu(x)$', 'Interpreter', 'LaTeX')

print('-depsc', '-tiff', '-r600', 'Normalised_computation_times_vs_input_arguments.eps')
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print('-dpdf', 'Normalised_computation_times_vs_input_arguments.pdf')
