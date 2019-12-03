%% Reducded-Order Modeling of Heat Conduction Problems

% Introduction:
% This project applying ROM on heat conduction problems. The original PDE
% has been shown in the SolveState.m. The result of ROM solution will be
% evaluated by comparing it with the solution from chebfun/pde15s.
%
% Heat Diffusion, Heat Advection-Deffusion and Burgers Probelm will be
% considered in this project.
%

%% Initialization
clear all
close all
clc
LW = 'LineWidth';

%% ======================= Part1: POD Computation =========================
Tf = 2;
t = 0:0.01:Tf;
x = chebfun(@(x) x,[-1,1]);
% Initial Condition
u0 = exp(-25*(x+0.5)^2);
% Parameters of the PDE
alpha = [0.1 0.1 0.1];
c = [0 1 0];
beta =[0 0 1];
% ==================== Compute POD modes ======================
% Case1: Heat Diffusion Problem(Linear)
u1 = SolveState(t,u0,alpha(1),beta(1),c(1));
[U1,S1,V1] = svd(u1);
for i=1:3
    plot(U1(:,i),LW,2);hold on
end
title('POD,Case1')
legend('mode1','mode2','mode3')

% Case2: Heat Advection-Diffusion Probelm(Linear)
u2 = SolveState(t,u0,alpha(2),beta(2),c(2));
[U2,S2,V2] = svd(u2);
figure
for i=1:3
    plot(U2(:,i),LW,2);hold on
end
title('POD,Case2')
legend('mode1','mode2','mode3')

% Case3: Burgers Probelm(No Linear)
u3 = SolveState(t,u0,alpha(3),beta(3),c(3));
[U3,S3,V3] = svd(u3);
figure
for i=1:3
    plot(U3(:,i),LW,2);hold on
end
title('POD,Case3')
legend('mode1','mode2','mode3')

% Plot singular values
delta1=zeros(1,25);
delta2=zeros(1,25);
delta3=zeros(1,25);
for i=1:25
    delta1(i)=S1(i,i);
end
for i=1:25
    delta2(i)=S2(i,i);
end
for i=1:25
    delta3(i)=S3(i,i);
end
figure
semilogy(delta1,LW,2);hold on
semilogy(delta2,LW,2)
semilogy(delta3,LW,2)
title('POD,Singular Values')
legend('Diffusion','Advection-Diffusion','Burgers')

% ===================== Least POD modes ==========================
%----case1----
% Frobenius Matrix Norm is applied to evaluate the energy contrained in
% each singular values
FN_u1 = sum(sum(u1.^2));
POD_u1 = 0;
for i=1:size(u1,2)
    POD_u1 = POD_u1+U1(:,i)*S1(i,i)*V1(:,i)';
    if sum(sum(POD_u1.^2))/FN_u1 > 0.999
        Least_mode1 = i;
        break
    end
end
fprintf(["Least Modes Required for 0.999 accuracy for case1:Diffusion is: %f\n"], Least_mode1)

%----case2----
FN_u2 = sum(sum(u2.^2));
POD_u2 = 0;
for i=1:size(u2,2)
    POD_u2 = POD_u2+U2(:,i)*S2(i,i)*V2(:,i)';
    if sum(sum(POD_u2.^2))/FN_u2 > 0.999
        Least_mode2 = i;
        break
    end
end
fprintf(["Least Modes Required for 0.999 accuracy for case2:Advection-Diffusion is: %f\n"], Least_mode2)

%----case3----
FN_u3 = sum(sum(u3.^2));
POD_u3 = 0;
for i=1:size(u3,2)
    POD_u3 = POD_u3+U3(:,i)*S3(i,i)*V3(:,i)';
    if sum(sum(POD_u3.^2))/FN_u3 > 0.999
        Least_mode3 = i;
        break
    end
end
fprintf(["Least Modes Required for 0.999 accuracy for case3:Burgers is: %f\n"], Least_mode3)

%% ========================== Part2: ROM =================================
% In this part, using ROM.m to transfer original PDEs into ODEs, and solve
% this ROM probelm with ODE solver. Then, compare the snapshot of ROM
% results with those corresponding snapshot from original PDE solutions

% Case1
figure;plot(u1(:,51),LW,2);hold on
r = 3;
dydt1 = @(t,Y) ROM(t,Y,U1,r,alpha(1),beta(1),c(1),x);
[t,Y] = ode23(dydt1,t,U1(:,1:r)'*u0);
urom=U1(:,1:r)*Y';
plot(urom(:,51),LW,2)

% Case2
figure;plot(u2(:,51),LW,2);hold on
r = 3;
dydt2 = @(t,Y) ROM(t,Y,U2,r,alpha(2),beta(2),c(2),x);
[t,Y] = ode23(dydt2,t,U2(:,1:r)'*u0);
urom=U2(:,1:r)*Y';
plot(urom(:,51),LW,2)

% Case3
figure;plot(u3(:,151),LW,2);hold on
r=3;
dydt3 = @(t,Y) ROM(t,Y,U3,r,alpha(3),beta(3),c(3),x);
[t,Y] = ode23(dydt3,t,U3(:,1:r)'*u0);
urom=U3(:,1:r)*Y';
plot(urom(:,151),LW,2)