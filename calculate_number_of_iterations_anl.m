function N = calculate_number_of_iterations_anl(j_T, n_T, m, method)

N = 0;

switch method
    
    case 'RFM'
        
        N = (1/6) * (j_T+1)*(j_T+2)*(j_T+3);
        
    case 'CEM-I'
        
        N = (1/6) * (j_T+1)*(j_T+2)*(j_T+3)*(n_T+1);
        
    case 'CEM-II'
        
        if j_T <= n_T
            N_phi = (1/48) * j_T*(j_T+1)*(j_T+2)*(j_T+5);
        else
            N_phi = (n_T/48) * (  j_T*(4*j_T^2+24*j_T+34) ...
                                - (n_T-1)*(n_T-2)*(n_T-5) ...
                                - 2*j_T*n_T*(3*j_T-2*n_T+12)  );
        end
        
        N = calculate_number_of_iterations_anl(j_T, n_T, m, 'CEM-I') - N_phi;
        
    case 'CEM-III'
        
        N_rho_1 = (1/24) * (j_T+1)*(j_T+2)*(j_T+3)*(n_T+1)*(3*j_T+2*n_T);
        
        if j_T <= n_T
            N_rho_2 = (1/48) * j_T*(j_T+1)*(j_T+2)*(j_T^2+5*j_T+1);
        else
            N_rho_2 = (n_T/48) * (  (n_T+1)*(n_T+2)*(n_T^2+5*n_T+1) ...
                                  + j_T*(3*j_T^3+16*j_T^2+20*j_T+6) ...
                                  - 2*n_T*(n_T^3+17*n_T+3) ...
                                  - j_T*n_T*(2*j_T^2+2*j_T*n_T-3*n_T^2+16*n_T-14)  );
        end
        
        N_rho = N_rho_1 - N_rho_2;
        
        N = (m-1)*calculate_number_of_iterations_anl(j_T, n_T, m, 'CEM-II') + N_rho;
        
    otherwise
        error('Incorrect method')
        
end
