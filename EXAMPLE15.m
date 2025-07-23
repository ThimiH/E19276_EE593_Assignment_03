%% Implementation of LMS based Beamforming

%% We will start with simulating a Colored Gaussian signal and use it in a
%% channel modeling problem implemented using LMS and varients of LMS

clear all
close all
clc

%% Simulation parameters

SIGMAa = 0.1 ; % Standard deviation of the desired signal
SIGMAb = 1.0 ; % Standard deviation of the jammer signal
OMEGA0 = 1.0 ;

TH0 = pi ;
PHI0 = pi*sin(TH0) ;

L = 10000 ; % Signal length

MU = 0.005 ;

%% Backup variables
Xbkp = zeros(2,L) ;
Dbkp = zeros(1,L) ;
Ybkp = zeros(1,L) ;
Ebkp = zeros(1,L) ;
Wbkp = zeros(2,L) ;

%% Initialization

n = 0 ;

w = [0 ; 0] ;


%% Run the system

for n = 1:L


    a = SIGMAa*randn(1) ;
    b = SIGMAb*randn(1) ;

    THa = 2*pi*rand(1) - pi ;
    THb = 2*pi*rand(1) - pi ;

    s = a*cos(OMEGA0*n + THa) ;
    v = b*cos(OMEGA0*n + THb) ;

    d = a*cos(OMEGA0*n + THa) + b*cos(OMEGA0*n + THb - PHI0) ;
    xn = a*cos(OMEGA0*n + THa) + b*cos(OMEGA0*n + THb) ;
    xt = a*sin(OMEGA0*n + THa) + b*sin(OMEGA0*n + THb) ;

    X = [xn ; xt] ;

    y = w'*X ;

    e = d - y ;

    w = w + 2*MU*e*X ; % LMS update equation

    %% Prepare for the next iteration

    % Backup variables
    Xbkp(:,n) = X ; % Update Xbkp
    Wbkp(:,n) = w ; % Update Wbkp
    Ybkp(:,n) = y ; % Update Ybkp
    Ebkp(:,n) = e ; % Update Ebkp
    Dbkp(:,n) = d ; % Update Dbkp

end

mean(e.^2)/mean(xn.^2)

figure(1) ; plot(1:L,Xbkp(1,:),1:L,Xbkp(2,:)) ; xlabel('n') ; ylabel('x(n)')
figure(2) ; plot(1:L,Dbkp) ; xlabel('n') ; ylabel('d(n)')
figure(3) ; plot(1:L,Ybkp) ; xlabel('n') ; ylabel('y(n)')
figure(4) ; plot(1:L,Ebkp) ; xlabel('n') ; ylabel('e(n)')
figure(5) ; plot(1:L,Wbkp(1,:),1:L,Wbkp(2,:))
figure(6)
th = -pi:(pi/100):pi ;
ph = pi*sin(th) ;
Gth = (cos(ph) - w(1)).^2 + (sin(ph) - w(2)).^2 ;
polar(th,Gth)