%% Implementation of LMS based Equalization

%% We will start with simulating a Colored Gaussian signal and use it in a
%% channel modeling problem implemented using LMS and varients of LMS

clear all
close all
clc

randn('seed',0)

%% Simulation parameters
% Channel model is given by G(z) = 2 / (1 + 0.8z^(-1) - 0.04z^(-2))
C = [1 ; -0.01 ; 0.5] ;
K = 0.5 ;

C2 = [1 ; -1.5 ; 0.54] ;
K2 = 0.6 ;

L = 3.0e5 ; % Signal length

Nsd = length(C) - 1 ; % System denominator order
Nsn = length(K) ; % System numerator order
Nw = 3 ; % Wiener Filter order

MU = 0.001 ;

%% Plot "Non-Adaptive" Equalized Channel Characteristics

CHAN = tf([K 0 0],C',1) 
CHAN2 = tf([K2 0 0],C2',1) 

ZFEQ = tf(C',[K 0 0],1)  

[H,W] = freqz([K 0 0],C') ;

figure(1)
subplot(2,1,1)
plot(W/pi,20*log10(abs(H)))
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
title('Channel Spectrum')
grid on
hold on

subplot(2,1,2)
plot(W/pi,atan2(imag(H),real(H))*180/pi)
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase (degrees)')
grid on
hold on

[H,W] = freqz([K2 0 0],C2') ;

figure(1)
subplot(2,1,1)
plot(W/pi,20*log10(abs(H)),'r')

subplot(2,1,2)
plot(W/pi,atan2(imag(H),real(H))*180/pi,'r')

figure(2)
freqz(C',[K 0 0])
title('Zero Forced Equalizer Spectrum')

ZFEQCHAN = series(CHAN,ZFEQ) 
[ZFEQCAHN_N,ZFEQCHAN_D] = tfdata(ZFEQCHAN,'v') ;

[H,W] = freqz(ZFEQCAHN_N,ZFEQCHAN_D) ;

figure(10)
subplot(2,1,1)
plot(W/pi,20*log10(abs(H)))
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Magnitude (dB)')
title('Equalizer Performance - Channel with Equalizer')
grid on
hold on

subplot(2,1,2)
plot(W/pi,atan2(imag(H),real(H))*180/pi)
xlabel('Normalized Frequency (\times\pi rad/sample)')
ylabel('Phase (degrees)')
grid on
hold on

ZFEQCHAN2 = series(CHAN2,ZFEQ) 
[ZFEQCAHN2_N,ZFEQCHAN2_D] = tfdata(ZFEQCHAN2,'v') ;

[H,W] = freqz(ZFEQCAHN2_N,ZFEQCHAN2_D) ;

figure(10)
subplot(2,1,1)
plot(W/pi,20*log10(abs(H)),'c')

subplot(2,1,2)
plot(W/pi,atan2(imag(H),real(H))*180/pi,'c')

% return

%% Backup variables
Vbkp = zeros(1,L) ;
Xbkp = zeros(1,L) ;
Dbkp = zeros(1,L) ;
Ybkp = zeros(1,L) ;
Ebkp = zeros(1,L) ;
Wbkp = zeros(Nw,L) ;

%% Initialization
v = 0 ;
V = v*eye(Nsn,1) ;
w = zeros(Nw,1) ;

Wbkp(:,1) = w ;

x = 0 ;
X = x*eye(Nw,1) ;
Xs = x*eye(Nsd,1) ;
Xbkp(:,1) = x ;

y = w'*X ;

d = v ;
e = d - y ;


%% Run the system "in real time"

T_C12 = 0.8e5 ; % Time when the channel changes
T_S12 = 0.5e5 ; % Time when first signal changes
T_S23 = 2.0e5 ; % Time when second signal changes

SD = 0.5 ; % 0.1 ; % 
SX = 0.01 ; % 0.1 ; % 

for I = 2:L

    y = w'*X ; % Wiener filter output
    d = v + 0.0*randn(1) ; % "Desired" signal (channel input plus noise)
    e = d - y ; % "Error" signal

    w = w + 2*MU*e*X ; % LMS update equation

    %% Prepare for the next iteration

    % Change Frequency Appropriately

    if (I<T_S12)
        FREQ = 1.5 ; 
    elseif (I<T_S23)
        FREQ = 30 ;
    else
        FREQ = 1.5 ;
    end
    
    v = 1.0*sin(2*pi*I*FREQ/1e4)+ SX*randn(1) ; % System input update equation
    
    for J = Nsn:-1:2
        V(J) = V(J-1) ;
    end
    V(1) = v ;
    
    % Change the channel
    if (I == T_C12)
    C = C2 ;
    K = K2 ;
    end
    
    x = ( K'*V - C(2:end)'*Xs)/C(1) ; % Update equation of the channel
    x = x + SX*randn(1) ; % Add noise to the channel
    
    % Update X (i.e. x -> X(1) -> X(2) -> X(3) etc. )
    for J = Nw:-1:2
        X(J) = X(J-1) ;
    end
    X(1) = x ;

    for J = Nsd:-1:2
        Xs(J) = Xs(J-1) ;
    end
    Xs(1) = x ;

    % Backup variables
    Xbkp(:,I) = x ; % Update Xbkp
    Wbkp(:,I) = w ; % Update Wbkp
    Ybkp(:,I) = y ; % Update Ybkp
    Ebkp(:,I) = e ; % Update Ebkp
    Dbkp(:,I) = d ; % Update Dbkp
    Vbkp(:,I) = v ; % Update Vbkp

end

ZFEQ_Y = lsim(ZFEQ,Xbkp) ;

figure(3)
plot(ZFEQ_Y)
xlabel('n') ; 
ylabel('y_{zero forced}(n)')

LMSEQ_N = Wbkp(:,T_C12-1)' ; 
LMSEQ_D = eye(size(LMSEQ_N)) ; % Something like [1 0 0 ... 0]
LMSEQ = tf(LMSEQ_N,LMSEQ_D,1) 

LMSEQCHAN = series(CHAN,LMSEQ) 
[LMSEQCAHN_N,LMSEQCHAN_D] = tfdata(LMSEQCHAN,'v') ;

[H,W] = freqz(LMSEQCAHN_N,LMSEQCHAN_D) ;

figure(10)
subplot(2,1,1)
plot(W/pi,20*log10(abs(H)),'r')

subplot(2,1,2)
plot(W/pi,atan2(imag(H),real(H))*180/pi,'r')

LMSEQ2_N = Wbkp(:,end)' ; 
LMSEQ2_D = eye(size(LMSEQ2_N)) ; % Something like [1 0 0 ... 0]
LMSEQ2 = tf(LMSEQ2_N,LMSEQ2_D,1) 

LMSEQCHAN2 = series(CHAN2,LMSEQ2) 
[LMSEQCAHN2_N,LMSEQCHAN2_D] = tfdata(LMSEQCHAN2,'v') ;

[H,W] = freqz(LMSEQCAHN2_N,LMSEQCHAN2_D) ;

figure(10)
subplot(2,1,1)
plot(W/pi,20*log10(abs(H)),'m')
hold off

subplot(2,1,2)
plot(W/pi,atan2(imag(H),real(H))*180/pi,'m')
hold off

figure(11) ; plot(1:L,Xbkp) ; xlabel('n') ; ylabel('x(n)')
title('Wiener Filter Input')
figure(12) ; plot(1:L,Dbkp) ; xlabel('n') ; ylabel('d(n)')
title('Wiener Filter Reference Input')
axis([0 L -max(abs(Dbkp)*1.5) max(abs(Dbkp)*1.5)])
figure(13) ; plot(1:L,Ybkp) ; xlabel('n') ; ylabel('y(n)')
title('Wiener Filter Output')
axis([0 L -max(abs(Dbkp)*1.5) max(abs(Dbkp)*1.5)])
figure(14) ; plot(1:L,Ebkp) ; xlabel('n') ; ylabel('e(n)')
title('Error Signal')
axis([0 L -2 2])
figure(15) ; plot(1:L,Wbkp(1,:),1:L,Wbkp(2,:),1:L,Wbkp(3,:))
xlabel('n') ; ylabel('W(n)')
legend('w_0(n)','w_1(n)','w_2(n)')
title('Tap Weights')

