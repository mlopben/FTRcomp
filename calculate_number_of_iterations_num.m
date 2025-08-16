function N = calculate_number_of_iterations_num(j_T, n_T, m, method)

N = 0;

switch method
    
    case 'RFM'
        
        for j = 0:j_T
            for k = 0:j
                for l = 0:k
                    N = N + 1;
                end
            end
        end
        
    case 'CEM-I'
        
        for j = 0:j_T
            for k = 0:j
                for l = 0:k
                    for n = 0:n_T
                        N = N + 1;
                    end
                end
            end
        end
        
    case 'CEM-II'
        
        for j = 0:j_T
            for k = 0:j
                for l = 0:k
                    n_0 = min(max(0,k-2*l),n_T);
                    for n = n_0:n_T
                        N = N + 1;
                    end
                end
            end
        end
        
    case 'CEM-III'
        
        for j = 0:j_T
            for k = 0:j
                for l = 0:k
                    n_0 = min(max(0,k-2*l),n_T);
                    for n = n_0:n_T
                        for p = 1:(j+m+n-1)
                            N = N + 1;
                        end
                    end
                end
            end
        end
        
    otherwise
        error('Incorrect method')
        
end
