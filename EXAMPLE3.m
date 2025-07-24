%% Steapest Descent Method: A two dimensional problem

%% TASK 1: By direct differentiation, find the local minima. (First show 
%% that R is positive definite by computing the eigenvalues of R). 

%% TASK 2: Change the step size parameter (Mu) and observe how it affects
%% the convergence. Obtain the range of Mu so that convergence happens.
%% Compare the maximum and minimum values of Mu with the eqigenvalues of R
%% (actually the reciprocal (1/a) of the eigenvalues). 

%% TASK 3: What happens if the initial conditions are changed? Set the
%% initial conditions of to be precisely on the eigenvectors of R and
%% observe what happens. Observe the convergence trajectory if the system
%% is initialized away fron the eigenvectors.

close all
clear all
clc

% Define the objective function
R = [4 -2 ; -1 3] ;
p = [2 ; 1] ;
D = 6 ;
f = @(x,y) R(1,1)*x.^2 + (R(1,2) + R(2,1))*x.*y  + R(2,2)*y.^2 - 2*p(1,1)*x - 2*p(2,1)*y + D; 

% Define the derivative of the objective function
Df = @(x,y) [(2*R(1,1)*x + (R(1,2) + R(2,1))*y - 2*p(1,1)) ; (2*R(2,2)*y + (R(1,2) + R(2,1))*x - 2*p(2,1))] ;

MaxIter = 1000 ;

EPS = 1e-3 ;

x = [3 ; 3] ; % Initial condition
xp = [inf ; inf] ;

X = x ;

Mu = 0.19 ; % Step size parameter

I = 0 ;

while ((norm(x-xp)>EPS)&(I<MaxIter))
I = I + 1 ;
xp = x ;
x = x - Mu*Df(x(1),x(2)) ;
X(:,end+1) = x ;
end

t = -5:0.1:5 ;
s = -4:0.1:4 ;
[T,S] = meshgrid(t,s) ;
Z = f(T,S) ;

figure(1)

subplot(3,2,[1 3])
mesh(t,s,Z)
hold on
plot3(X(1,:),X(2,:),f(X(1,:),X(2,:)),'.r-')
xlabel('x')
ylabel('y')
zlabel('f(x,y)')
axis vis3d
hold off

subplot(3,2,5)
contourf(t,s,Z)
hold on
plot(X(1,:),X(2,:),'.r-')
xlabel('x')
ylabel('y')
axis equal
axis tight
hold off

subplot(3,2,2)
plot(0:(size(X,2)-1),X(1,:),'.-')
xlabel('Iteration Number (k)')
ylabel('x_k')
grid on

subplot(3,2,4)
plot(0:(size(X,2)-1),X(2,:),'.-')
xlabel('Iteration Number (k)')
ylabel('y_k')
grid on

subplot(3,2,6)
plot(0:(size(X,2)-1),f(X(1,:),X(2,:)),'.-')
xlabel('Iteration Number (k)')
ylabel('f(x_k,y_k)')
grid on
