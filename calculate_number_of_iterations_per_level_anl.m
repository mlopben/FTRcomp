function N_iter = calculate_number_of_iterations_per_level_anl(method, level, j_T, n_T, m)

switch level
    
    case 'j'
        N_iter = j_T+1;
    
    case 'k'
        N_iter = (1/2) * (j_T+1)*(j_T+2);
    
    case 'l'
        N_iter = calculate_number_of_iterations_anl(j_T, n_T, m, 'RFM');
    
    case 'n'
        switch method
            case 'RFM'
                error('RFM does not have this level of iteration')
            case 'CEM-I'
                N_iter = calculate_number_of_iterations_anl(j_T, n_T, m, 'CEM-I');
            case 'CEM-II'
                N_iter = calculate_number_of_iterations_anl(j_T, n_T, m, 'CEM-II');
            case 'CEM-III'
                N_iter = calculate_number_of_iterations_anl(j_T, n_T, m, 'CEM-II');
            otherwise
                error('Wrong method [A]')
        end
    
    case 'p'
        switch method
            case 'CEM-III'
                N_iter = calculate_number_of_iterations_anl(j_T, n_T, m, 'CEM-III');
            otherwise
                error('Wrong method [B]')
        end
    
    otherwise
        error('Incorrect level of iteration')
end
