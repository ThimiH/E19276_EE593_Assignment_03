%% Steapest Descent Method: A one dimensional problem

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

x = 3 ; % Initial condition
xp = -inf ;

X = x ;

MU = 0.9; % Step size parameter

I = 0 ;

while ((abs(x-xp)>EPS)&(I<MaxIter))
I = I + 1 ;
xp = x ;
x = x - MU*Df(x) ;
X(end+1) = x ;
end

t = -5:0.01:4 ;
y = f(t) ;

figure(1)

subplot(2,2,[1 3])
plot(t,y,'b-',X,f(X),'r.-')
axis([-5 4 0 25])
xlabel('x')
ylabel('f(x)')

subplot(2,2,2)
plot(0:(length(X)-1),X,'.-')
xlabel('Iteration Number (k)')
ylabel('x_k')

subplot(2,2,4)
plot(0:(length(X)-1),f(X),'.-')
xlabel('Iteration Number (k)')
ylabel('f(x_k)')
