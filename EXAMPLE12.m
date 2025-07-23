%% "Real time" Implementation of VLMS

%% We will start with simulating a Colored Gaussian signal and use it in a
%% channel modeling problem implemented using LMS and varients of LMS

clear all
close all
clc

%% Simulation parameters
% Channel model is given by G(z) = 2 - 1.2z^(-1) + 0.8z^(-2)
C = [2 ; -1.2 ; 0.8] ;

Sig = 0.5 ; % Standard deviation of the Wiener process
L = 2000 ; % Signal length

Ns = length(C) ; % System order
Nw = 3 ; % Wiener Filter order

MU = zeros(Nw,1) ; 
mu_max = 0.2 ; % Maximum step
mu_min = 1e-3 ; % Minimum step
RHO = 0.01 ; % Step size of step update

%% Backup variables
Xbkp = zeros(1,L) ;
Dbkp = zeros(1,L) ;
Ybkp = zeros(1,L) ;
Ebkp = zeros(1,L) ;
Wbkp = zeros(Nw,L) ;
Mbkp = zeros(Nw,L) ;

%% Initialization
x = 0.1 ;
w = zeros(Nw,1) ;

X = x*eye(Nw,1) ;
Xp = 0*ones(size(X)) ;
Xs = x*eye(Ns,1) ;
Xbkp(:,1) = x ;
Wbkp(:,1) = w ;
Mbkp(:,1) = MU ;

d = C'*Xs ;
y = w'*X ;
e = d - y ; 
ep = 0 ; 

%% Run the system

for I = 2:L

    d = C'*Xs ; % "Desired" signal (channel output)
    
    y = w'*X ; % Wiener filter output
    
    e = d - y ; % "Error" signal
    
    g = e*X ; 
    g_1 = ep*Xp ;
    
    MU = MU + RHO*g.*g_1 ; % Update step size
    %MU = MU + RHO*sign(g).*sign(g_1) ; % Update step size
    
    MU = min(max(MU,mu_min),mu_max) ; % Limit step size
    
    w = w + 2*MU.*g ; % LMS update equation  

%% Prepare for the next iteration

    x = randn(1) ; % System input update equation
    
    % Update X (i.e. x -> X(1) -> X(2) -> X(3) etc. )
    for J = Nw:-1:2
        X(J) = X(J-1) ;
    end
    X(1) = x ;

    % Update Xs
    for J = Ns:-1:2
        Xs(J) = Xs(J-1) ;
    end
    Xs(1) = x ;

    % Backup variables
    Xbkp(:,I) = x ; % Update Xbkp
    Wbkp(:,I) = w ; % Update Wbkp
    Ybkp(:,I) = y ; % Update Wbkp
    Ebkp(:,I) = e ; % Update Wbkp
    Dbkp(:,I) = d ; % Update Wbkp
    Mbkp(:,I) = MU ; % Update Wbkp
    ep = e ;
    Xp = X ;

end


figure(1) ; plot(1:L,Xbkp) ; xlabel('n') ; ylabel('x(n)')
figure(2) ; plot(1:L,Dbkp) ; xlabel('n') ; ylabel('d(n)')
figure(3) ; plot(1:L,Ybkp) ; xlabel('n') ; ylabel('y(n)')
figure(4) ; plot(1:L,Ebkp) ; xlabel('n') ; ylabel('e(n)')
figure(5) ; plot(1:L,Wbkp(1,:),1:L,Wbkp(2,:),1:L,Wbkp(3,:))
xlabel('n') ; ylabel('W(n)')
legend('w_1(n)','w_2(n)','w_3(n)')
figure(6) ; plot(1:L,Mbkp(1,:),1:L,Mbkp(2,:),1:L,Mbkp(3,:))
xlabel('n') ; ylabel('\mu(n)')
legend('\mu_1(n)','\mu_2(n)','\mu_3(n)')