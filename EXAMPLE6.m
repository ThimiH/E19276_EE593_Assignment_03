%% An Application Using Simulated Signals

%% TASK 1: Change the bandwidth parameter 'a' and observe the effect.

%% TASK 2: For each 'a' change the 'Mu' and observe the effect.

%% TASK 3: For each 'a' and 'Mu' change the initialization of 'w' and
%% observe the effect. 


close all
clear all
clc

%% Parameters

a = 0.5 ;

L = 100000 ;

df = 2/L ;
f = 0:df:(2 - df) ;

%% Generate signals

v = wgn(1,L,0) ; 

v = v - mean(v) ;
v = v/std(v) ;

G = tf(sqrt([1-a^2 0]),[1 -a],1) 

x = lsim(G,v) ;
x = x + 0.0*randn(size(x)) ;

P = tf([1 -4],[1 0],1)

d = lsim(P,x) ;
d = d + 0.0*randn(size(d)) ;

%% Plot signal spectra

figure(1)

subplot(4,1,1)
plot(f,abs(fft(v)).^2 /L)
xlabel('Notmalized Frequency (\times \pi rad)')
ylabel('\Phi_{vv}')

subplot(4,1,2)
plot(f,abs(fft(x)).^2 /L)
xlabel('Notmalized Frequency (\times \pi rad)')
ylabel('\Phi_{xx}')

subplot(4,1,3)
plot(f,abs(fft(d)).^2 /L)
xlabel('Notmalized Frequency (\times \pi rad)')
ylabel('\Phi_{dd}')

subplot(4,1,4)
plot(f,abs(fft(d).*conj(fft(x))) /L)
xlabel('Notmalized Frequency (\times \pi rad)')
ylabel('\Phi_{dx}')

%% Calculate signal statistics

x1 = x(1:(end-1))' ;
x0 = x(2:end)' ;

R = [sum(x0.^2) sum(x0.*x1) ; sum(x1.*x0) sum(x1.^2)] /L 

d1 = d(1:(end-1))' ;
d0 = d(2:end)' ;

p = [sum(x0.*d0) ; sum(x1.*d0)]/(L-1)

%% Use LMS to calculate wo

w = [4;4] ;
w_bkp = w ;
wp = [inf;inf] ;

EPS = 1e-6 ;

Mu = 0.01 ;

MaxIter = 1000 ;

I = 0 ;

W = w ;

while(I<MaxIter)%((norm(w-wp)>EPS)&(I<MaxIter))
    I = I + 1 ;
    wp = w ;
    y = w'*[x0(I) ; x1(I)] ;
    e = d0(I) - y ;
    w = w + 2*Mu*e*[x0(I) ; x1(I)] ;
    W(:,end+1) = w ;
end

wo_WH = [1 a ; a 1]\[1-4*a ; a-4] ;
disp('Ideal solution from Wiener-Hopf Equation')
disp(wo_WH)

disp('Calculated solution')
disp(w)

%% Output

y = w' * [x0 ; x1] ;

e = d0 - y ;

[eVec,eVal] = eig(R) ;

% Sort the eigenvectors of R in the descending order
[eVal,Idx] = sort(diag(eVal),'descend') ;
eVec = eVec(:,Idx) 
eVal = eVal

%% Output signal spectra

figure(2)
subplot(3,1,1)
plot(f(1:(end-1))*(L-1)/L,abs(fft(y)).^2 /L)
xlabel('Notmalized Frequency (\times \pi rad)')
ylabel('\Phi_{yy}')

subplot(3,1,2)
plot(f(1:(end-1))*(L-1)/L,abs(fft(e).*conj(fft(y))) /L)
xlabel('Notmalized Frequency (\times \pi rad)')
ylabel('\Phi_{ey}')

subplot(3,1,3)
plot(f(1:(end-1))*(L-1)/L,abs(fft(e)).^2 /L)
xlabel('Notmalized Frequency (\times \pi rad)')
ylabel('\Phi_{ee}')

%% Performance surfaces

ww0 = -10:0.1:10 ;
ww1 = -15:0.1:5 ;
[WW0,WW1] = meshgrid(ww0,ww1) ;

XIXI = sum(d0.^2)/(L-1) - 2*WW0*p(1) - 2*WW1*p(2) ...
    + WW0.*WW0*R(1,1) + WW0.*WW1.*(R(1,2)+R(2,1)) ...
    + WW1.*WW1*R(2,2) ;

XI = sum(d0.^2)/(L-1) - 2*W(1,:)*p(1) - 2*W(2,:)*p(2) ...
    + W(1,:).*W(1,:)*R(1,1) + W(1,:).*W(2,:).*(R(1,2)+R(2,1)) ...
    + W(2,:).*W(2,:)*R(2,2) ; 

figure(3) ; 

mesh(ww0,ww1,XIXI)
hold on
plot3(W(1,:),W(2,:),XI,'.r-')
xlabel('w_0')
ylabel('w_1')
zlabel('\xi_k')
title('performance surface \xi(w)')
axis vis3d
hold off

figure(4) ; 
contourf(ww0,ww1,XIXI)
hold on
plot(W(1,:),W(2,:),'.r-')
xlabel('w_0')
ylabel('w_1')
title('Level Curves of \xi(w)')
% The following 7 lineas are for plotting the eigenvectors on the level
% curves
T = 5 ;
EV1p = wo_WH + T*eVec(:,1) ;
EV1n = wo_WH - T*eVec(:,1) ;
EV2p = wo_WH + T*eVec(:,2) ;
EV2n = wo_WH - T*eVec(:,2) ;
plot([EV1p(1) EV1n(1)],[EV1p(2) EV1n(2)],'g-')
plot([EV2p(1) EV2n(1)],[EV2p(2) EV2n(2)],'m-')
text(EV1n(1),EV1n(2),['\lambda =' num2str(eVal(1))],'Color','green')
text(EV2n(1),EV2n(2),['\lambda =' num2str(eVal(2))],'Color','magenta')
axis equal
axis tight
hold off

figure(5) ; 
subplot(2,3,1)
plot(0:(size(W,2)-1),W(1,:),'.-')
xlabel('Iteration Number (k)')
ylabel('w_0(k)')
grid on

subplot(2,3,2)
plot(0:(size(W,2)-1),W(1,:)*eVec(1,1)+W(2,:)*eVec(2,1),'.-')
xlabel('Iteration Number (k)')
ylabel(['v' char(39) '{}_1(k)'])
grid on

subplot(2,3,3)
plot(0:(size(W,2)-1),W(2,:),'.-')
xlabel('Iteration Number (k)')
ylabel('w_1(k)')
grid on

subplot(2,3,4)
plot(0:(size(W,2)-1),W(1,:)*eVec(1,2)+W(2,:)*eVec(2,2),'.-')
xlabel('Iteration Number (k)')
ylabel(['v' char(39) '{}_2(k)'])
grid on

subplot(2,3,5)
plot(0:(size(W,2)-1),XI,'.-')
xlabel('Iteration Number (k)')
ylabel('\xi_k')
grid on

subplot(2,3,6)
plot(0:(size(W,2)-1),sqrt(W(1,:).^2 + W(2,:).^2),'.-')
xlabel('Iteration Number (k)')
ylabel('|| w_k||')
grid on
