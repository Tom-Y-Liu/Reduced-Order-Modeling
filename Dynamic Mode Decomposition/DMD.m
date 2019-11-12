%% Reduced-Order Modeling - Dynamic Mode Decomposition

% Introduction:
% This project appling Dynamic Mode Decomposition(DMD) to an 'unknown'
% dynamic system. And do forecast based on the DMD system.
%

%% Initialization
clear all
clc
close all
LW = 'LineWidth';
load('Data.mat')
u = D;
t = 0:0.001:4;
dt = t(2) - t(1);

%% ===================== Part1: Compute the DMD ==========================
r = 21;% dimension of the DMD mode
X1 = u(:,1:2000);
X2 = u(:,2:2001);
[U,S,V] = svd(X1,'econ');
Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);

% A is n*n. Too large.
A = X2*pinv(X1);

% Atilde is r*r
Atilde = Ur'*X2*Vr*inv(Sr);

% Compute the components of the solution of DMD
[W,D] = eig(Atilde);
phi=X2*Vr*inv(Sr)*W;
lambda = diag(D);
omega = log(lambda)/dt;

% rank the DMD modes
[Y,I] = sort(real(omega),'descend');
lambda = lambda(I);
omega = omega(I);
phi = phi(:,I);

% ==== Back to high dimensional space ====
% build u_DMD
b=phi\X1(:,1);
time_dynamics=zeros(r,length(t));
for iter=1:length(t)
    time_dynamics(:,iter) = (b.*exp(omega*t(iter)));
end
u_DMD = phi*time_dynamics;

%% ==================== Part2: Plot DMD modes ============================
for i = [2 4 5]
    plot(real(phi(:,i)),LW,2);hold on
end
for i = [1 3]
    plot(real(phi(:,i)),'--',LW,2);hold on
end
legend('mode2','mode4','mode5','mode1','mode3')
title('DMD modes (r=5)')

%% ===================== Part3: Plot eigenvalues =========================
figure
for i = 1:r
    plot(real(lambda(i)),imag(lambda(i)),'x');hold on
end
title('eigenvalues')
xlabel('Real')
ylabel('Imag')

%% ========================= Part4: Forecast =============================
figure;
plot(u(:,4001),LW,2);hold on
plot(real(u_DMD(:,4001)),'--',LW,2)
title('t=4.0')
legend('Truth','DMD forecast')
