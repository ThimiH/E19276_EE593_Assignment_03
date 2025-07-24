%% Analysis of Steepest Descent Method from EXAMPLE3.m
%% Testing different step sizes and their effect on convergence

close all;
clear all;
clc;

%% Problem Setup (from EXAMPLE3.m)
% Define the objective function
R = [4 -2; -1 3];
p = [2; 1];
D = 6;
f = @(x,y) R(1,1)*x.^2 + (R(1,2) + R(2,1))*x.*y + R(2,2)*y.^2 - 2*p(1,1)*x - 2*p(2,1)*y + D;

% Define the derivative of the objective function
Df = @(x,y) [(2*R(1,1)*x + (R(1,2) + R(2,1))*y - 2*p(1,1)); 
             (2*R(2,2)*y + (R(1,2) + R(2,1))*x - 2*p(2,1))];

%% Analysis of R matrix
eigenvalues = eig(R);
fprintf('Eigenvalues of R: %.4f, %.4f\n', eigenvalues(1), eigenvalues(2));
fprintf('R is positive definite: %s\n', string(all(eigenvalues > 0)));

[eigenvectors, eigenvals_matrix] = eig(R);
fprintf('Eigenvectors of R (columns):\n');
disp(eigenvectors);

% Optimal solution
w_star = R \ p;
fprintf('Optimal w (local minima): [%.4f, %.4f]\n', w_star(1), w_star(2));
fprintf('Minimum function value: %.4f\n\n', f(w_star(1), w_star(2)));

%% Parameters
MaxIter = 1000;
EPS = 1e-3;
initial_condition = [2; -1]; % Initial condition

% Different step sizes to test
step_sizes = [-0.01, 0, 0.01, 0.1, 0.19, 0.20, 0.21];

% Grid for plotting
t = -5:0.1:5;
s = -4:0.1:4;
[T, S] = meshgrid(t, s);
Z = f(T, S);

%% Analyze each step size
results = struct();

for idx = 1:length(step_sizes)
    Mu = step_sizes(idx);
    fprintf('Testing step size Mu = %.3f\n', Mu);
    
    % Initialize
    x = initial_condition;
    xp = [inf; inf];
    X = x;
    I = 0;
    converged = false;
    
    % Steepest descent iterations
    while ((norm(x - xp) > EPS) && (I < MaxIter))
        I = I + 1;
        xp = x;
        
        if Mu == 0
            % Special case: no update when step size is 0
            break;
        else
            x = x - Mu * Df(x(1), x(2));
        end
        
        X(:, end+1) = x;
        
        % Check for divergence
        if any(~isfinite(x)) || norm(x) > 1e6
            fprintf('  Diverged at iteration %d\n', I);
            break;
        end
    end
    
    % Check convergence
    if norm(x - xp) <= EPS
        converged = true;
        fprintf('  Converged in %d iterations\n', I);
    else
        fprintf('  Did not converge (reached max iterations or diverged)\n');
    end
    
    % Store results
    results(idx).Mu = Mu;
    results(idx).iterations = I;
    results(idx).converged = converged;
    results(idx).trajectory = X;
    results(idx).final_point = x;
    results(idx).final_error = norm(x - w_star);
    
    %% Create and save figure
    figure('Position', [100, 100, 432, 252]); % 6" x 3.5" at 72 DPI
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [12 7]);
    set(gcf, 'PaperPosition', [0 0 12 7]);
    
    % Main contour plot with trajectory
    subplot(2, 3, [1 4]);
    contourf(t, s, Z, 20);
    hold on;
    plot(X(1, :), X(2, :), 'r.-', 'LineWidth', 2, 'MarkerSize', 8);
    plot(w_star(1), w_star(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    plot(initial_condition(1), initial_condition(2), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    colorbar;
    xlabel('x');
    ylabel('y');
    title(sprintf('Steepest Descent: \\mu = %.3f', Mu));
    axis equal;
    axis tight;
    grid on;
    
    % Add text with iteration count
    if converged
        text(0.02, 0.98, sprintf('Converged in %d iterations', I), ...
             'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'white', 'EdgeColor', 'black');
    else
        text(0.02, 0.98, sprintf('Did not converge (%d iterations)', I), ...
             'Units', 'normalized', 'VerticalAlignment', 'top', ...
             'BackgroundColor', 'white', 'EdgeColor', 'black');
    end
    hold off;
    
    % x trajectory
    subplot(2, 3, 2);
    plot(0:(size(X, 2)-1), X(1, :), '.-');
    xlabel('Iteration Number (k)');
    ylabel('x_k');
    title('x component');
    grid on;
    
    % y trajectory
    subplot(2, 3, 3);
    plot(0:(size(X, 2)-1), X(2, :), '.-');
    xlabel('Iteration Number (k)');
    ylabel('y_k');
    title('y component');
    grid on;
    
    % Function value trajectory
    subplot(2, 3, [5 6]);
    f_values = f(X(1, :), X(2, :));
    plot(0:(size(X, 2)-1), f_values, '.-');
    xlabel('Iteration Number (k)');
    ylabel('f(x_k, y_k)');
    title('Function Value');
    grid on;
    
    % Save figure
    mu_str = sprintf('%.2f', abs(Mu));
    if Mu < 0
        mu_str = ['neg' mu_str];
    elseif Mu == 0
        mu_str = '0.00';
    end
    mu_str = strrep(mu_str, '.', '_');
    
    filename = sprintf('Answer_1_c2_Fig_MU_%s.png', mu_str);
    print(filename, '-dpng', '-r300');
    fprintf('  Saved: %s\n\n', filename);
    
    close(gcf);
end

%% Summary Analysis
fprintf('=== CONVERGENCE SUMMARY ===\n');
fprintf('Step Size    | Iterations | Converged | Final Error\n');
fprintf('-------------|------------|-----------|------------\n');

for idx = 1:length(results)
    fprintf('%9.3f    |%10d  |%9s  |%10.4f\n', ...
            results(idx).Mu, results(idx).iterations, ...
            string(results(idx).converged), results(idx).final_error);
end

%% Theoretical Analysis
fprintf('\n=== THEORETICAL ANALYSIS ===\n');
fprintf('For convergence, step size must satisfy: 0 < μ < 2/λ_max\n');
fprintf('Maximum eigenvalue λ_max = %.4f\n', max(eigenvalues));
fprintf('Theoretical upper bound: μ < %.4f\n', 2/max(eigenvalues));
fprintf('Minimum eigenvalue λ_min = %.4f\n', min(eigenvalues));
fprintf('For fastest convergence: μ ≈ 2/(λ_max + λ_min) = %.4f\n', ...
        2/(max(eigenvalues) + min(eigenvalues)));

fprintf('\nAnalysis complete! Figures saved with naming convention: Answer_1_c_Fig_MU_xx.png\n');
