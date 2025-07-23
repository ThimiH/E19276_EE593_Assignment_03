%% "Real time" Implementation of Signed-Regressor LMS

%% We will start with simulating a Colored Gaussian signal and use it in a
%% channel modeling problem implemented using LMS and varients of LMS

clear all
close all
clc

%% Simulation parameters
% Channel model is given by G(z) = 2 - 1.2z^(-1) + 0.8z^(-2)
C = [2 ; -1.2 ; 0.8] ;

Sig = 0.5 ; % Standard deviation of the Wiener process
L = 200 ; % Signal length

Ns = length(C) ; % System order
Nw = 3 ; % Wiener Filter order

MU = 0.1 ; 
PSI = 0.01 ;

%% Backup variables
Xbkp = zeros(1,L) ;
Dbkp = zeros(1,L) ;
Ybkp = zeros(1,L) ;
Ebkp = zeros(1,L) ;
Wbkp = zeros(Nw,L) ;

%% Initialization
x = 0.1 ;
w = zeros(Nw,1) ;

X = x*eye(Nw,1) ;
Xs = x*eye(Ns,1) ;
Xbkp(:,1) = x ;
Wbkp(:,1) = w ;

d = C'*Xs ;
y = w'*X ;
e = d - y ; 

%% Run the system

for I = 2:L

    d = C'*Xs ; % "Desired" signal (channel output)
    
    y = w'*X ; % Wiener filter output
    
    e = d - y ; % "Error" signal
    
    w = w + 2*MU*e*sign(X) ; % LMS update equation  

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

end

figure(1)
subplot(5,1,1) ; plot(1:L,Xbkp) ; xlabel('n') ; ylabel('x(n)')
subplot(5,1,2) ; plot(1:L,Dbkp) ; xlabel('n') ; ylabel('d(n)')
subplot(5,1,3) ; plot(1:L,Ybkp) ; xlabel('n') ; ylabel('y(n)')
subplot(5,1,4) ; plot(1:L,Ebkp) ; xlabel('n') ; ylabel('e(n)')
subplot(5,1,5) ; plot(1:L,Wbkp(1,:),1:L,Wbkp(2,:),1:L,Wbkp(3,:))
xlabel('n') ; ylabel('W(n)')
legend('w_1(n)','w_2(n)','w_3(n)')