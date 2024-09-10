% MATLAB script for computing E[ξN] and E[(ξN - E[ξN])^2]

% Run the initial script
run('synthetic_data.m'); %
disp('Now running main_script...');

% Parameters
A0_values = Largest_Oms_As;   % amplitude, ensure this is defined in 'synthetic_data.m'
omegas = [omega_gK, omega_gK1]; % multiple omega values, ensure these are defined
sigma = 1; % Standard deviation of random variable epsilon
Delta = 0.1; % Time step interval

% Set the seed for reproducibility
set_seed = 123;
rng(set_seed);

% Array of different N values (number of time steps) to test
N_values = [100, 500, 1000, 10000]; % N values

% Pre-allocate arrays for results
E_xiN_all_N = zeros(length(N_values) * length(A0_values), length(omegas)); % Expected values for each N and amplitude
Var_xiN_all_N = zeros(length(N_values) * length(A0_values), length(omegas)); % Variance for each N and amplitude
result_counter = 1;

% Store amplitude information
amplitudes = [];

% Loop over different N values
for n_idx = 1:length(N_values)
    N = N_values(n_idx); % Set the current N value
    
    % Time vector for the current N
    t = Delta * (1:N);
    
    % Loop over both amplitude values
    for j = 1:length(A0_values)
        A0 = A0_values(j);
        disp(['Running with amplitude A0 = ', num2str(A0), ' and N = ', num2str(N)]);
        
        % Pre-allocate arrays for the current N and amplitude
        E_xiN = zeros(1, length(omegas));
        Var_xiN = zeros(1, length(omegas));
        
        % Loop over omega values
        for i = 1:length(omegas)
            omega = omegas(i);
            fprintf('Omega %d: %f, Amplitude: %f, N = %d\n', i, omega, A0, N);
            
            % Generate epsilon for this iteration
            epsilon = sigma * randn(1, N); % Generate noise with standard normal distribution
            
            % Compute ξN for the current N and amplitude
            xiN = (1/N) * sum((A0 * sin(omega * t)) + epsilon);
            
            % Store expected value and variance for the current N and amplitude
            E_xiN(i) = xiN; % Accumulate expected value
            Var_xiN(i) = (xiN - mean(E_xiN)).^2; % Accumulate variance
        end
        
        % Store the results for this N and amplitude
        E_xiN_all_N(result_counter, :) = E_xiN;
        Var_xiN_all_N(result_counter, :) = Var_xiN;
        amplitudes(result_counter, 1) = A0; % Store the amplitude used
        result_counter = result_counter + 1;
    end
end

% Combine N values and amplitudes into the result table
E_xiN_all_N_2 = [repelem(N_values', length(A0_values)), amplitudes, E_xiN_all_N];
Var_xiN_all_N_2 = [repelem(N_values', length(A0_values)), amplitudes, Var_xiN_all_N];

% Create tables with amplitude included
Table_exp_mean = array2table(E_xiN_all_N_2, 'VariableNames', {'Sample_size', 'Amplitude', 'E[ξN] with Omega_1', 'E[ξN] with Omega_2'});
Table_variance = array2table(Var_xiN_all_N_2, 'VariableNames', {'Sample_size', 'Amplitude', 'V[ξN] with Omega_1', 'V[ξN] with Omega_2'});

% Display the tables
fprintf('Expected Values of ξN for each omega, N, and amplitude:\n');
fprintf('------------------------------------\n');
disp(Table_exp_mean)

fprintf('Variance of ξN for each omega, N, and amplitude:\n');
fprintf('------------------------------------\n');
disp(Table_variance)
