
clear
clc

% j_T = 10;
j_T = 40;
n_T = 15;
m = 3;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  NUMBER OF ITERATIONS OF RFM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
N_sim = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            N_sim = N_sim + 1;
        end
    end
end

% Analytic
N_anl = (1/6) * (j_T + 1) * (j_T + 2) * (j_T + 3);

% Symbolic calculation
syms sym_j_T j k l
N_sym = symsum(symsum(symsum(1, l, 0, k), k, 0, j), j, 0, sym_j_T);
sym_j_T = j_T;
N_sym = double(subs(N_sym));

fprintf('\n\nNumber of iterations for RFM:')
fprintf('\n-------------------------------------------------')
fprintf('\n   Sim    : %g', N_sim);
fprintf('\n   Anl    : %g', N_anl);
fprintf('\n   Sym    : %g', N_sym);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  NUMBER OF ITERATIONS OF CEM-I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation
N_sim = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            for n = 0:n_T
                N_sim = N_sim + 1;
            end
        end
    end
end

% Analytic
N_anl = (1/6) * (j_T + 1) * (j_T + 2) * (j_T + 3) * (n_T + 1);

% Symbolic calculation
syms sym_j_T sym_n_T j k l n
N_sym = symsum(symsum(symsum(symsum(1, n, 0, sym_n_T), l, 0, k), k, 0, j), j, 0, sym_j_T);
sym_j_T = j_T;
sym_n_T = n_T;
N_sym = double(subs(N_sym));

fprintf('\n\nNumber of iterations for CEM-I:')
fprintf('\n-------------------------------------------------')
fprintf('\n   Sim    : %g', N_sim);
fprintf('\n   Anl    : %g', N_anl);
fprintf('\n   Sym    : %g', N_sym);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  NUMBER OF ITERATIONS OF CEM-II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation A
N_sim_A = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            n_0 = min(max(0,k-2*l),n_T);
            for n = n_0:n_T
                N_sim_A = N_sim_A + 1;
            end
        end
    end
end

% Simulation B
N_sim_B = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            n_0 = min(max(0,k-2*l),n_T);
            N_sim_B = N_sim_B + n_T - n_0 + 1;
        end
    end
end

% Simulation C
N_sim_C1 = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            n_0 = min(max(0,k-2*l),n_T);
            N_sim_C1 = N_sim_C1 + n_T + 1;
        end
    end
end
N_sim_C2 = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            n_0 = min(max(0,k-2*l),n_T);
            N_sim_C2 = N_sim_C2 + n_0;
        end
    end
end
N_sim_C = N_sim_C1 - N_sim_C2;

% Simulation D & E
N_sim_D = (1/6) * (j_T + 1) * (j_T + 2) * (j_T + 3) * (n_T + 1); % This is the number of iterations of CEM-I = N_Sim_C1
N_sim_E = (1/6) * (j_T + 1) * (j_T + 2) * (j_T + 3) * (n_T + 1); % This is the number of iterations of CEM-I = N_Sim_C1
D = 0;
E = 0;
for j = 0:j_T
    for k = 0:j
        % fprintf('\n j = %d, k = %d -----------------', j, k)
        for l = 0:k
            
            n_0 = min(max(0,k-2*l),n_T);
            D = D + n_0;
            
            if l <= floor((k-n_T)/2)
                % fprintf('\nlow: l = %d', l)
                E = E + n_T;
            end
            if l >= (floor((k-n_T)/2)+1) && l <= floor(k/2)
                % fprintf('\nmid: l = %d', l)
                E = E + k - 2*l;
            end
            
        end
    end
end
N_sim_D = N_sim_D - D;
N_sim_E = N_sim_E - E;

% Simulation F
N_sim_F = (1/6) * (j_T + 1) * (j_T + 2) * (j_T + 3) * (n_T + 1); % This is the number of iterations of CEM-I = N_Sim_C1
X = 0;
for j = 0:j_T
    for k = 0:j
        
        if k <= n_T
            for l = 0 : floor(k/2)
                X = X + k - 2*l;
            end
            for l = (floor(k/2)+1) : k
                X = X + 0;
            end
        end
        
        if k >= n_T + 1
            for l = 0 : floor((k-n_T)/2)
                X = X + n_T;
            end
            for l = (floor((k-n_T)/2)+1) : floor(k/2)
                X = X + k - 2*l;
            end
            for l = (floor(k/2)+1) : k
                X = X + 0;
            end
        end
        
    end
end
N_sim_F = N_sim_F - X;

% Simulation G
N_sim_G = (1/6) * (j_T + 1) * (j_T + 2) * (j_T + 3) * (n_T + 1); % This is the number of iterations of CEM-I = N_Sim_C1
X = 0;
if j_T <= n_T % ---------------------- This is for the case j_T <= n_T
    % fprintf('\nj_T <= n_T\n')
    for j = 0:j_T
        for k = 0:j
            for l = 0:floor(k/2)
                X = X + k - 2*l;
            end
        end
    end
else % ------------------------------- This is for the case j_T > n_T
    % fprintf('\nj_T > n_T\n')

    for j = 0:n_T
        for k = 0:j
            for l = 0:floor(k/2)
                X = X + k - 2*l;
            end
        end
    end
    
    for j = (n_T+1):j_T
        for k = 0:n_T
            for l = 0:floor(k/2)
                X = X + k - 2*l;
            end
        end
        for k = (n_T+1):j
            for l = 0:floor((k-n_T)/2)
                X = X + n_T;
            end
            for l = (floor((k-n_T)/2)+1):floor(k/2)
                X = X + k - 2*l;
            end
        end
    end
end
N_sim_G = N_sim_G - X;

% Analytic
if j_T <= n_T % ---------------------- This is for the case j_T <= n_T
    X = ( (1/48) * j_T * (j_T + 1) * (j_T + 2) * (j_T + 5) );
else % ------------------------------- This is for the case j_T > n_T
    X = ( (n_T/48) * (  j_T*(4*j_T^2+24*j_T+34)  -  (n_T-1)*(n_T-2)*(n_T-5)  -  2*j_T*n_T*(3*j_T-2*n_T+12) ) );
end
N_anl = ( (1/6) * (j_T + 1) * (j_T + 2) * (j_T + 3) * (n_T + 1) )   -   X;

% Symbolic calculation
syms sym_j_T sym_n_T j k l n
if j_T <= n_T % ---------------------- This is for the case j_T <= n_T
    X = symsum(symsum(symsum(k-2*l, l, 0, k/2), k, 0, j), j, 0, sym_j_T);
else % ------------------------------- This is for the case j_T > n_T
    X = symsum(symsum(symsum(k-2*l, l, 0, k/2), k, 0, j), j, 0, sym_n_T) ...
      + symsum(symsum(symsum(k-2*l, l, 0, k/2), k, 0, sym_n_T) + symsum(symsum(sym_n_T, l, 0, ((k-sym_n_T)/2)) + symsum(k-2*l, l, ((k-sym_n_T)/2)+1, k/2), k, sym_n_T+1, j), j, sym_n_T+1, sym_j_T);
end
sym_j_T = j_T;
sym_n_T = n_T;
N_sym = (1/6) * (j_T + 1) * (j_T + 2) * (j_T + 3) * (n_T + 1) - double(subs(X));

fprintf('\n\nNumber of iterations for CEM-II ')
if j_T <= n_T
    fprintf('(j_T <= n_T) :')
else
    fprintf('(j_T > n_T) :')
end
fprintf('\n-------------------------------------------------')
fprintf('\n   Sim A  : %g', N_sim_A);
fprintf('\n   Sim B  : %g', N_sim_B);
fprintf('\n   Sim C  : %g', N_sim_C);
fprintf('\n   Sim D  : %g', N_sim_D);
fprintf('\n   Sim E  : %g', N_sim_E);
fprintf('\n   Sim F  : %g', N_sim_F);
fprintf('\n   Sim G  : %g', N_sim_G);
fprintf('\n   Anl    : %g', N_anl);
fprintf('\n   Sym    : %g', N_sym);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  NUMBER OF ITERATIONS OF CEM-III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulation A
N_sim_A = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            n_0 = min(max(0,k-2*l),n_T);
            for n = n_0:n_T
                for p = 1:(j+m+n-1)
                    N_sim_A = N_sim_A + 1;
                end
            end
        end
    end
end

% Simulation B
N_CEM_2 = 0;
Y = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            n_0 = min(max(0,k-2*l),n_T);
            for n = n_0:n_T
                N_CEM_2 = N_CEM_2 + 1;
                Y = Y + j + n;
            end
        end
    end
end
N_sim_B = (m-1)*N_CEM_2 + Y;

% Simulation C
Y1 = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            for n = 0:n_T
                Y1 = Y1 + j + n;
            end
        end
    end
end
Y2 = 0;
for j = 0:j_T
    for k = 0:j
        for l = 0:k
            n_0 = min(max(0,k-2*l),n_T);
            Y2 = Y2 + n_0*(n_0+2*j-1)/2;
        end
    end
end
N_sim_C = (m-1)*N_CEM_2 + Y1 - Y2;

% Analytic
Y1 = ((3*j_T + 2*n_T)*(j_T + 1)*(j_T + 2)*(j_T + 3)*(n_T + 1))/24;
if j_T <= n_T % ---------------------- This is for the case j_T <= n_T
    Y2 = (j_T*(j_T + 1)*(j_T + 2)*(j_T^2 + 5*j_T + 1))/48;
    X = ( (1/48) * j_T * (j_T + 1) * (j_T + 2) * (j_T + 5) );
else % ------------------------------- This is for the case j_T > n_T
    Y2 = (n_T/48) * (  (n_T + 1)*(n_T + 2)*(n_T^2 + 5*n_T + 1) ...
                     +  j_T*(3*j_T^3 + 16*j_T^2 + 20*j_T + 6) ...
                     -  2*n_T*(n_T^3 + 17*n_T + 3) ...
                     -  j_T*n_T*(2*j_T^2 + 2*j_T*n_T - 3*n_T^2 + 16*n_T - 14)  );
    X = ( (n_T/48) * (  j_T*(4*j_T^2+24*j_T+34)  -  (n_T-1)*(n_T-2)*(n_T-5)  -  2*j_T*n_T*(3*j_T-2*n_T+12) ) );
end
N_CEM_2 = ( (1/6) * (j_T + 1) * (j_T + 2) * (j_T + 3) * (n_T + 1) )   -   X;
N_anl = (m-1)*N_CEM_2 + Y1 - Y2;

% Symbolic calculation
syms sym_j_T sym_n_T j k l n
Y1 = symsum(symsum(symsum(symsum(j+n, n, 0, sym_n_T), l, 0, k), k, 0, j), j, 0, sym_j_T);
if j_T <= n_T % ---------------------- This is for the case j_T <= n_T
	Y2 = symsum(symsum(symsum((k-2*l)*(k-2*l+2*j-1)/2, l, 0, k/2), k, 0, j), j, 0, sym_j_T);
else % ------------------------------- This is for the case j_T > n_T
    Y2 = symsum(symsum(symsum((k-2*l)*(k-2*l+2*j-1)/2, l, 0, k/2), k, 0, j), j, 0, sym_n_T) ...
       + symsum( ...
                 symsum(symsum((k-2*l)*(k-2*l+2*j-1)/2, l, 0, k/2), k, 0, sym_n_T) ...
               + symsum(symsum(sym_n_T*(sym_n_T+2*j-1)/2, l, 0, ((k-sym_n_T)/2)) + symsum((k-2*l)*(k-2*l+2*j-1)/2, l, ((k-sym_n_T)/2)+1, k/2), k, sym_n_T+1, j), j, sym_n_T+1, sym_j_T);
end
sym_j_T = j_T;
sym_n_T = n_T;
Y = Y1 - Y2;
N_sym = (m-1)*N_CEM_2 + double(subs(Y));

fprintf('\n\nNumber of iterations for CEM-III ')
if j_T <= n_T
    fprintf('(j_T <= n_T) :')
else
    fprintf('(j_T > n_T) :')
end
fprintf('\n-------------------------------------------------')
fprintf('\n   Sim A  : %g', N_sim_A);
fprintf('\n   Sim B  : %g', N_sim_B);
fprintf('\n   Sim C  : %g', N_sim_C);
fprintf('\n   Anl    : %g', N_anl);
fprintf('\n   Sym    : %g', N_sym);

fprintf('\n\n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  PLOT THE NUMBER OF ITERATIONS FOR THE CEM MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

j_T_values = 0:60;
n_T_values = 0:30;

method = 'CEM-I';
N_iter = zeros(numel(j_T_values), numel(n_T_values));
for j_T = j_T_values
    for n_T = n_T_values
        N_iter(j_T+1, n_T+1) = calculate_number_of_iterations_num(j_T, n_T, -1, method);
    end
end
figure
surf(n_T_values, j_T_values, N_iter)
xlabel('n_T')
ylabel('j_T')
zlabel('N_{iter}')
title(method)

method = 'CEM-II';
N_iter = zeros(numel(j_T_values), numel(n_T_values));
for j_T = j_T_values
    for n_T = n_T_values
        N_iter(j_T+1, n_T+1) = calculate_number_of_iterations_num(j_T, n_T, -1, method);
    end
end
figure
surf(n_T_values, j_T_values, N_iter)
xlabel('n_T')
ylabel('j_T')
zlabel('N_{iter}')
title(method)

method = 'CEM-III';
m = 1;
N_iter = zeros(numel(j_T_values), numel(n_T_values));
for j_T = j_T_values
    for n_T = n_T_values
        N_iter(j_T+1, n_T+1) = calculate_number_of_iterations_num(j_T, n_T, m, method);
    end
end
figure
surf(n_T_values, j_T_values, N_iter)
xlabel('n_T')
ylabel('j_T')
zlabel('N_{iter}')
title(sprintf('%s (m=%d)', method, m))

method = 'CEM-III';
m = 2;
N_iter = zeros(numel(j_T_values), numel(n_T_values));
for j_T = j_T_values
    for n_T = n_T_values
        N_iter(j_T+1, n_T+1) = calculate_number_of_iterations_num(j_T, n_T, m, method);
    end
end
figure
surf(n_T_values, j_T_values, N_iter)
xlabel('n_T')
ylabel('j_T')
zlabel('N_{iter}')
title(sprintf('%s (m=%d)', method, m))

method = 'CEM-III';
m = 3;
N_iter = zeros(numel(j_T_values), numel(n_T_values));
for j_T = j_T_values
    for n_T = n_T_values
        N_iter(j_T+1, n_T+1) = calculate_number_of_iterations_num(j_T, n_T, m, method);
    end
end
figure
surf(n_T_values, j_T_values, N_iter)
xlabel('n_T')
ylabel('j_T')
zlabel('N_{iter}')
title(sprintf('%s (m=%d)', method, m))
