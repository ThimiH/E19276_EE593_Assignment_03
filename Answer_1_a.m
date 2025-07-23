%% Analysis of Fixed Point Iteration Methods from EXAMPLE1.m
%% Testing convergence behavior for different initial values

clear all
close all
clc

%% Parameters
MAX_ITER = 1000;
EPS = 1e-3;

% Range of initial values to test
x0_range = -2:0.1:4;  % Test from -2 to 4 with step 0.1

% Initialize arrays to store results
iterations_method1 = zeros(size(x0_range));
iterations_method2 = zeros(size(x0_range));
iterations_method3 = zeros(size(x0_range));
iterations_method4 = zeros(size(x0_range));

converged_method1 = zeros(size(x0_range));
converged_method2 = zeros(size(x0_range));
converged_method3 = zeros(size(x0_range));
converged_method4 = zeros(size(x0_range));

fprintf('Analyzing convergence for %d initial values...\n', length(x0_range));

%% Test each initial value
for idx = 1:length(x0_range)
    x0 = x0_range(idx);
    
    %% Method 1: x = (x^2 + 1)/2
    try
        xp = -inf;
        x = x0;
        I = 0;
        
        while ((abs(x - xp) > EPS) && (I < MAX_ITER))
            I = I + 1;
            xp = x;
            x = (x^2 + 1)/2;
            
            % Check for NaN or Inf
            if ~isfinite(x)
                break;
            end
        end
        
        iterations_method1(idx) = I;
        converged_method1(idx) = (abs(x - xp) <= EPS) && isfinite(x);
    catch
        iterations_method1(idx) = MAX_ITER;
        converged_method1(idx) = 0;
    end
    
    %% Method 2: x = x^2 - x + 1
    try
        xp = -inf;
        x = x0;
        I = 0;
        
        while ((abs(x - xp) > EPS) && (I < MAX_ITER))
            I = I + 1;
            xp = x;
            x = x^2 - x + 1;
            
            % Check for NaN or Inf
            if ~isfinite(x)
                break;
            end
        end
        
        iterations_method2(idx) = I;
        converged_method2(idx) = (abs(x - xp) <= EPS) && isfinite(x);
    catch
        iterations_method2(idx) = MAX_ITER;
        converged_method2(idx) = 0;
    end
    
    %% Method 3: x = sqrt(2x - 1)
    try
        xp = -inf;
        x = x0;
        I = 0;
        
        while ((abs(x - xp) > EPS) && (I < MAX_ITER))
            I = I + 1;
            xp = x;
            
            % Check if argument is valid for real square root
            if (2*x - 1) >= 0
                x = sqrt(2*x - 1);
            else
                x = sqrt(2*x - 1); % Will produce complex result
            end
            
            % Check for NaN or Inf
            if ~isfinite(x)
                break;
            end
        end
        
        iterations_method3(idx) = I;
        converged_method3(idx) = (abs(x - xp) <= EPS) && isfinite(x);
    catch
        iterations_method3(idx) = MAX_ITER;
        converged_method3(idx) = 0;
    end
    
    %% Method 4: x = 2 - (1/x)
    try
        xp = -inf;
        x = x0;
        I = 0;
        
        % Skip if x0 is zero (division by zero)
        if x0 == 0
            iterations_method4(idx) = MAX_ITER;
            converged_method4(idx) = 0;
        else
            while ((abs(x - xp) > EPS) && (I < MAX_ITER))
                I = I + 1;
                xp = x;
                
                % Check for division by zero
                if x == 0
                    break;
                end
                
                x = 2 - (1/x);
                
                % Check for NaN or Inf
                if ~isfinite(x)
                    break;
                end
            end
            
            iterations_method4(idx) = I;
            converged_method4(idx) = (abs(x - xp) <= EPS) && isfinite(x);
        end
    catch
        iterations_method4(idx) = MAX_ITER;
        converged_method4(idx) = 0;
    end
end

%% Display Results Summary
fprintf('\n=== CONVERGENCE ANALYSIS SUMMARY ===\n');
fprintf('Initial value range: [%.1f, %.1f] with step %.1f\n', min(x0_range), max(x0_range), x0_range(2)-x0_range(1));
fprintf('Maximum iterations allowed: %d\n', MAX_ITER);
fprintf('Convergence tolerance: %.0e\n\n', EPS);

methods = {'Method 1: x = (x^2 + 1)/2', 'Method 2: x = x^2 - x + 1', ...
           'Method 3: x = sqrt(2x - 1)', 'Method 4: x = 2 - (1/x)'};
converged = {converged_method1, converged_method2, converged_method3, converged_method4};
iterations = {iterations_method1, iterations_method2, iterations_method3, iterations_method4};

for m = 1:4
    fprintf('%s:\n', methods{m});
    conv_count = sum(converged{m});
    fprintf('  Converged cases: %d/%d (%.1f%%)\n', conv_count, length(x0_range), 100*conv_count/length(x0_range));
    
    if conv_count > 0
        conv_idx = find(converged{m});
        fprintf('  Convergent range: x0 ∈ [%.1f, %.1f]\n', min(x0_range(conv_idx)), max(x0_range(conv_idx)));
        fprintf('  Average iterations (converged cases): %.1f\n', mean(iterations{m}(conv_idx)));
        fprintf('  Min/Max iterations: %d/%d\n', min(iterations{m}(conv_idx)), max(iterations{m}(conv_idx)));
    else
        fprintf('  No convergent cases found in the tested range\n');
    end
    fprintf('\n');
end

%% Plotting Results (Only Converging Cases)
figure(1);
hold on;

% Method 1 - only converging cases
conv_idx1 = find(converged_method1);
if ~isempty(conv_idx1)
    plot(x0_range(conv_idx1), iterations_method1(conv_idx1), 'r.', 'MarkerSize', 18, 'DisplayName', 'Method 1: (x^2+1)/2');
end

% Method 2 - only converging cases
conv_idx2 = find(converged_method2);
if ~isempty(conv_idx2)
    plot(x0_range(conv_idx2), iterations_method2(conv_idx2), 'b.', 'MarkerSize', 18, 'DisplayName', 'Method 2: x^2-x+1');
end

% Method 3 - only converging cases
conv_idx3 = find(converged_method3);
if ~isempty(conv_idx3)
    plot(x0_range(conv_idx3), iterations_method3(conv_idx3), 'm.', 'MarkerSize', 18, 'DisplayName', 'Method 3: sqrt(2x-1)');
end

% Method 4 - only converging cases
conv_idx4 = find(converged_method4);
if ~isempty(conv_idx4)
    plot(x0_range(conv_idx4), iterations_method4(conv_idx4), 'g.', 'MarkerSize', 18, 'DisplayName', 'Method 4: 2-1/x');
end

xlabel('Initial Value (x_0)');
ylabel('Iterations to Converge');
title('Convergence Analysis: Iterations vs Initial Value (Converging Cases Only)');
legend('Location', 'northwest');
grid on;

% Set y-limit based on actual converging iterations
all_conv_iters = [];
if ~isempty(conv_idx1), all_conv_iters = [all_conv_iters, iterations_method1(conv_idx1)]; end
if ~isempty(conv_idx2), all_conv_iters = [all_conv_iters, iterations_method2(conv_idx2)]; end
if ~isempty(conv_idx3), all_conv_iters = [all_conv_iters, iterations_method3(conv_idx3)]; end
if ~isempty(conv_idx4), all_conv_iters = [all_conv_iters, iterations_method4(conv_idx4)]; end

if ~isempty(all_conv_iters)
    ylim([0 max(all_conv_iters) + 5]);
else
    ylim([0 50]);
end
hold off;

%% Convergence Map
figure(2);
convergence_matrix = [converged_method1; converged_method2; converged_method3; converged_method4];
imagesc(x0_range, 1:4, convergence_matrix);
colormap([1 0 0; 0 1 0]); % Red for non-convergent, Green for convergent
colorbar('Ticks', [0.25, 0.75], 'TickLabels', {'No Convergence', 'Convergence'});
yticks(1:4);
yticklabels({'Method 1', 'Method 2', 'Method 3', 'Method 4'});
xlabel('Initial Value (x_0)');
ylabel('Method');
title('Convergence Map (Green = Converges, Red = Diverges)');
grid on;

%% Detailed Analysis for Specific Initial Values
fprintf('=== DETAILED ANALYSIS FOR SPECIFIC INITIAL VALUES ===\n');
test_values = [0.5, 0.9, 1.0, 1.1, 1.5, 2.0];

for x0_test = test_values
    fprintf('\nInitial value x0 = %.1f:\n', x0_test);
    
    % Find closest index in our range
    [~, idx] = min(abs(x0_range - x0_test));
    
    for m = 1:4
        if converged{m}(idx)
            fprintf('  %s: Converged in %d iterations\n', methods{m}, iterations{m}(idx));
        else
            fprintf('  %s: Did not converge (>%d iterations)\n', methods{m}, MAX_ITER);
        end
    end
end

fprintf('\nAnalysis complete!\n');