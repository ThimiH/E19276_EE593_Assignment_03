%% Steepest Descent Method: A one dimensional problem

%% TASK 1: Change the step size parameter (MU) and observe how it affects
%% the convergence. Obtain the range of Mu so that convergence happens.

%% TASK 2: What happens if the initial conditions are changed?

close all
clear all
clc

% Define the objective function
f = @(x) x.^2 + 2*x + 3 ; 

% Define the derivative of the objective function
Df = @(x) 2*x + 2 ;

MaxIter = 1000 ;
EPS = 1e-3 ;

% Range of MU values to test
MU_values = [1.1, 1, 0.9, 0.51, 0.5, 0.4, 0.1, 0.001, 0, -0.1];

% Test different MU values
fprintf('Testing Steepest Descent with different step sizes (MU values):\n');
fprintf('=============================================================\n');

for i = 1:length(MU_values)
    MU = MU_values(i);
    fprintf('\nTesting MU = %.1f:\n', MU);
    
    try
        runSteepestDescent(f, Df, MU, MaxIter, EPS, 3);
    catch ME
        fprintf('Error with MU = %.1f: %s\n', MU, ME.message);
    end
end

fprintf('\nAll figures saved successfully!\n');

% Function to run steepest descent and create plots
function runSteepestDescent(f, Df, MU, MaxIter, EPS, x0)
    x = x0; % Initial condition
    xp = -inf;
    X = x;
    I = 0;
    
    while ((abs(x-xp) > EPS) & (I < MaxIter))
        I = I + 1;
        xp = x;
        x = x - MU * Df(x);
        X(end+1) = x;
        
        % Check for divergence
        if abs(x) > 1000
            break;
        end
    end
    
    % Create the plot
    t = -5:0.01:4;
    y = f(t);
    
    % Set figure size to 6" width and 3.5" height
    figure('Units', 'inches', 'Position', [1, 1, 6, 3.5]);
    
    subplot(2,2,[1 3])
    plot(t, y, 'b-', 'LineWidth', 1.5)
    hold on
    plot(X, f(X), 'r.-', 'MarkerSize', 8, 'LineWidth', 1.2)
    hold off
    axis([-5 4 0 25])
    xlabel('x', 'FontSize', 10)
    ylabel('f(x)', 'FontSize', 10)
    title(sprintf('Objective Function (MU = %.1f)', MU), 'FontSize', 10)
    grid on
    legend('f(x) = x^2 + 2x + 3', 'Optimization Path', 'Location', 'best')
    
    subplot(2,2,2)
    plot(0:(length(X)-1), X, '.-', 'MarkerSize', 6, 'LineWidth', 1.2)
    xlabel('Iteration Number (k)', 'FontSize', 10)
    ylabel('x_k', 'FontSize', 10)
    title('Parameter Evolution', 'FontSize', 10)
    grid on
    
    subplot(2,2,4)
    plot(0:(length(X)-1), f(X), '.-', 'MarkerSize', 6, 'LineWidth', 1.2)
    xlabel('Iteration Number (k)', 'FontSize', 10)
    ylabel('f(x_k)', 'FontSize', 10)
    title('Function Value Evolution', 'FontSize', 10)
    grid on
    
    % Add overall title
    sgtitle(sprintf('Steepest Descent Method: MU = %.1f, Iterations = %d', MU, I), 'FontSize', 12)
    
    % Adjust layout
    tight_layout = true;
    if tight_layout
        % Adjust subplot spacing
        set(gcf, 'PaperUnits', 'inches');
        set(gcf, 'PaperSize', [15 5]);
        set(gcf, 'PaperPosition', [0 0 6 3.5]);
    end
    
    % Save the figure
    filename = sprintf('Answer_1_b_Fig_MU_%.1f.png', MU);
    saveas(gcf, filename, 'png');
    fprintf('Saved figure: %s (Iterations: %d, Final x: %.4f)\n', filename, I, X(end));
    
    % Return convergence info
    converged = (abs(x - xp) <= EPS) && (I < MaxIter);
    return
end
