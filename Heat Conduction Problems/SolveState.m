function u=SolveState(t,u0,alpha,beta,c)
% This function is from ROM course by Dr.Hessam Babaee of University of
% pittsburgh in 2019.

Domain = [-1 1];
% Make the right-hand side of the PDE.
pdefun = @(t,x,u) alpha*diff(u,2)-c*diff(u)- beta*u.*diff(u) + f(t,x);

% Assign boundary conditions.
bc.left = 'dirichlet';
bc.right = 'dirichlet';

% bc = 'periodic';
% Construct a chebfun of the space variable on the domain,
x = chebfun(@(x) x, Domain);
%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-6 , 'Ylim', [-3,3]);
%
%% Call pde15s to solve the problem.
[t, u] = pde15s(pdefun, t, u0, bc, opts);
