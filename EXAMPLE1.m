%% Solving x^2 - 2x + 1 = 0 using fixed point iteration

%% TASK: Change the initial condition x0 and observe the reult given by
%% each method. Find the range of initial conditions for which the
%% convergence happens. 

clear all
close all
clc

x0 = 0.9 ;

MAX_ITER = 1000 ;

%% Method 1: x = (x^2 + 1)/2

xp = -inf ;
x = x0 ;

X1 = x ;

EPS = 1e-3 ;

I = 0 ;

while ((abs(x - xp)>EPS) & (I<MAX_ITER))
    I = I+1 ;
    xp = x ;
    x = (x.^2 + 1)/2 ;
    X1(end+1) = x ;
end

%% Method 2: x = x^2 - x + 1

xp = -inf ;
x = x0 ;

X2 = x ;

EPS = 1e-3 ;

I = 0 ;

while ((abs(x - xp)>EPS) & (I<MAX_ITER))
    I = I+1 ;
    xp = x ;
    x = x.^2 - x + 1 ;
    X2(end+1) = x ;
end


%% Method 3: x = sqrt(1-2x)

xp = -inf ;
x = x0 ;

X3 = x ;

EPS = 1e-3 ;

I = 0 ;

while ((abs(x - xp)>EPS) & (I<MAX_ITER))
    I = I+1 ;
    xp = x ;
    x = sqrt(2*x -1) ;
    X3(end+1) = x ;
end

%% Method 4: x = 2 - (1/x)

xp = -inf ;
x = x0 ;

X4 = x ;

EPS = 1e-3 ;

I = 0 ;

while ((abs(x - xp)>EPS) & (I<MAX_ITER))
    I = I+1 ;
    xp = x ;
    x = 2 - (1./x) ;
    X4(end+1) = x ;
end

%% Plot the "convergence" of all methods
figure(1)
plot(0:(length(X1)-1),X1,'r.-')
xlabel('Iteration Number (k)')
ylabel('x_k')
grid on

figure(2)
plot(0:(length(X2)-1),X2,'b.-')
xlabel('Iteration Number (k)')
ylabel('x_k')
grid on 

figure(3)
plot3(real(X3),imag(X3),0:(length(X3)-1),'m.-')
xlabel('Real(x_k)')
ylabel('Imag(x_k)')
zlabel('Iteration Number (k)')
grid on 

figure(4)
plot(0:(length(X4)-1),X4,'b.-')
xlabel('Iteration Number (k)')
ylabel('x_k')
grid on