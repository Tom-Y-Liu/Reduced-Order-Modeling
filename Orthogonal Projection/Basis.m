function U=Basis(N,Identifier)
% BASIS can return 5 kind of bases. N is the dimension of the basis, and
% 'Identifier' is the corresponding name of each basis. 
%
% This function is from ROM course by Dr.Hessam Babaee of University of
% pittsburgh in 2019.
%

x = chebfun('x',[-1,1]);
X = linspace(-1,1,N);
U = repmat(x,1,N);
dx = 2/(N-1);
switch Identifier
    case 'Const'
        X = linspace(-1,1,N+1);
        for i=1:N
            U(:,i) = heaviside(x-X(i))-heaviside(x-X(i+1));
        end
    case 'Linear'
        for i=1:N
            U(:,i) =  max(1- abs(x-X(i))/dx,0);
        end
    case 'Legendre'
        for i=1:N
            U(:,i) = legpoly(i-1);
        end
    case 'Monomial'
        for i=1:N
            U(:,i) = x.^(i-1);
        end    
    case 'Fourier'
        U(:,1) = x-x+1;
         for i=2:2:N
            U(:,i)   = sin(i/2*pi*x);
            U(:,i+1) = cos(i/2*pi*x);
         end
        if (mod(N,2)==0)
            U = U(:,1:N);
        end
end