%% Implementation of LMS phase and amplitude estimation

%% We will start with simulating a Colored Gaussian signal and use it in a
%% channel modeling problem implemented using LMS and varients of LMS

clear all
close all
clc

%% Simulation parameters

OMEGA = 0.5 ;

a = 1 
TH0 = 3*pi/4 

L = 50000 ; % Signal length

MU = 0.1 ;

%% Backup variables
Xbkp = zeros(2,L) ;
Dbkp = zeros(1,L) ;
Ybkp = zeros(1,L) ;
Ebkp = zeros(1,L) ;
Wbkp = zeros(2,L) ;

%% Initialization

n = 0 ;
X = [sin(OMEGA*n) ; 0] ;
d = a*sin(OMEGA*n + TH0) + 0.01*randn(1) ;
w = [0 ; 0] ;


%% Run the system

for n = 1:L

    y = w'*X ;

    e = d - y ;

    w = w + 2*MU*e*X ; % LMS update equation

    %% Prepare for the next iteration

    X(2) = X(1) ; 
    X(1) = sin(OMEGA*n) ;
    
    d = a*sin(OMEGA*n + TH0) + 0.01*randn(1) ;
    
    % Backup variables
    Xbkp(:,n) = X ; % Update Xbkp
    Wbkp(:,n) = w ; % Update Wbkp
    Ybkp(:,n) = y ; % Update Ybkp
    Ebkp(:,n) = e ; % Update Ebkp
    Dbkp(:,n) = d ; % Update Dbkp

end

%% y(n) = w0*sin(OM*n) + w1*sin(OM*(n-1)) 
%%      = w0*sin(OM*n) + w1*(cos(OM))*(sin(OM*n)) - w1*(sin(OM))*(cos(OM*n))
%%      = (w0 + w1*cos(OM))sin(OM*n) - (w1*sin(OM))*cos(OM*n)
%%      = A(cos(PHI))sin(OM*n) + B(sin(PHI))cos(OM*n)

a_est = sqrt((w(1) + w(2)*cos(OMEGA))^2 + (-w(2)*sin(OMEGA))^2)
TH_est = atan2((-w(2)*sin(OMEGA)),(w(1) + w(2)*cos(OMEGA)))

figure(1) ; plot(1:L,Xbkp(1,:),1:L,Xbkp(2,:)) ; xlabel('n') ; ylabel('x(n)')
figure(2) ; plot(1:L,Dbkp) ; xlabel('n') ; ylabel('d(n)')
figure(3) ; plot(1:L,Ybkp) ; xlabel('n') ; ylabel('y(n)')
figure(4) ; plot(1:L,Ebkp) ; xlabel('n') ; ylabel('e(n)')
figure(5) ; plot(1:L,Wbkp(1,:),1:L,Wbkp(2,:))
