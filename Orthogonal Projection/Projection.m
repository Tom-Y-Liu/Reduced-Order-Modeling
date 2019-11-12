function pf=Projection(f,N,Identifier)
%PROJECTION returns the projection of function f onto 'Identifier' with N dimensional

U = Basis(N,Identifier);
M = U'*U;% The mass matrix
pf = U*(inv(M)*U'*f);