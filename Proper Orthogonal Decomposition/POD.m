%% Reduced-Order Modeling - Proper Othogonal Decomposition(POD)

% Introductions
% -------------
%
% This project is based on a time-dependent heat equation(PDE).
%
% Firstly,'chebfun' library will be applied to get the 'exact' solution of
% the PDE. Then, POD algorithm will be applied to transfer the PDE into a
% series of ODEs.
%
% Plots will be generated to compare the solution of POD and 'exact'
% solution.
%

%% Initialization
clc
clear all
close all

%% ======================== Part1:Solve the PDE ==========================
% This part solve the PDE with pde15s in chebfun.
%
% This section is from ROM course by Dr.Hessam Babaee of University of
% pittsburgh in 2019.
%

%% ------- Space-Time Domains -----------
Tf = 2;
t = 0:.01:Tf;
dom = [-1 1];
%% ------- Problem Parameters -----------
alpha = .1;

% Make the right-hand side of the PDE.
pdefun = @(t,x,u) alpha*diff(u,2);
% Assign boundary conditions.
bc.left = 'dirichlet';
bc.right = 'dirichlet';

x = chebfun(@(x) x, dom);
% and of the initial condition.
u0 = exp(-25*(x-.5).^2)+exp(-25*(x+.4).^2);
%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-6 , 'Ylim', [-3,3]);
%
%% Call pde15s to solve the problem.
[t, u] = pde15s(pdefun, t, u0, bc, opts);

waterfall(u)
title('u')


%% ============= Part2: Compute and Visualize the POD modes ==============
[U, S, V]=svd(u);
figure
plot(U(:,1));hold on
plot(U(:,2));
plot(U(:,3));
legend('Mode1','Mode2','Mode3')
title('First Three Modes')

%% ========== Part3: Plot the singular values of the first 25 modes ======
delta = zeros(25,1);
for i=1:25
    delta(i) = S(i,i);
end
figure
semilogy(delta);
title('Singular Values of First 25 modes')

%% ============ Part4: Check the number of modes needed ==================
FN_u = sum(sum(u.^2));
POD_u = 0;
for i=1:size(u,2)
    POD_u = POD_u+U(:,i)*S(i,i)*V(:,i)';
    if sum(sum(POD_u.^2))/FN_u > 0.999
        Least_mode = i;
        break
    end
end
fprintf(["Least Modes Required for 0.999 accuracy is: %f\n"], Least_mode)

%% ===================== Part5: Solve the ODEs ===========================
% This part solve the ODEs derived by POD process, and compare the
% solutions to the 'exact' solutions

LW = 'LineWidth';
figure
plot(u(:,101),LW,2);hold on
title('t=1')
for r=[2,3,4]
    Ur = U(:,1:r);
    M = Ur'*Ur;
    Y = zeros(r,size(t,1));
    %Initial Condition
    Y0 = exp(-25*(x-0.5)^2)+exp(-25*(x+0.4)^2);
    %d2U/dx2 secend order derivative
    diff2Ur=diff(Ur,2);
    %dY/dt = alpha*D*Y
    B = alpha*inv(M)'*(diff2Ur'*Ur)';
    %Solve the ROM
    [Q,lam]=eig(B);
    for j=1:size(t,1)
        Y(:,j) = Q*expm(t(j)*lam)*inv(Q)*(Ur'*Y0);
    end
    uROM = Ur*Y;
    plot(uROM(:,101),LW,2);
end
legend('PDE','ROM r=2','ROM r=3','ROM r=4')

figure
plot(u(:,201),LW,2);hold on
title('t=2')
for r=[2,3,4]
    Ur = U(:,1:r);
    M = Ur'*Ur;
    Y = zeros(r,size(t,1));
    %Initial Condition
    Y0 = exp(-25*(x-0.5)^2)+exp(-25*(x+0.4)^2);
    %d2U/dx2 secend order derivative
    diff2Ur=diff(Ur,2);
    %dY/dt = alpha*D*Y
    B = alpha*inv(M)'*(diff2Ur'*Ur)';
    %Solve the ROM
    [Q,lam]=eig(B);
    for j=1:size(t,1)
        Y(:,j) = Q*expm(t(j)*lam)*inv(Q)*(Ur'*Y0);
    end
    uROM = Ur*Y;
    plot(uROM(:,201),LW,2);
end
legend('PDE','ROM r=2','ROM r=3','ROM r=4')

figure
waterfall(uROM)
title('uROM r=4')