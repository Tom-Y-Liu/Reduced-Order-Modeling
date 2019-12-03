function g = ROM(t,Y,U,r,alpha,beta,c,x)
% ROM returns the right hand side of the ODE probelm, which is the
% projection of the original PDE probelm into basis Ur.

U=U(:,1:r);
% Mass Matrix
M=U'*U;
% Linear Part
L = -alpha*diff(U)'*diff(U)-c*U'*diff(U);
% Nonlinear Part
for i=1:r
    for j=1:r
        for k=1:r
            N(i,j,k) = -beta*(U(:,i).*diff(U(:,j)))'*U(:,k);
        end
    end
end
NY = zeros(r,1);
for k=1:r
    NY(k,1) = sum(N(:,:,k).*(Y*Y'),'all');
end
g=inv(M)*(NY+L*Y+U'*f(t,x));

%%
% I was tried to reduce the process of computing N to a O(n^2) probelm. But
% this method is not feasible to chefun based matrix, because U(:) works
% differently in matlab and chebfun.
%
% N = zeros(r,r,r);
% for k=1:r
%     A0 = [U(:,1)];
%     B = [];
%     C = [];
%     C0 = repmat(U(:,k),1,r);
%     for i=1:r
%         A0 = vertcat(A0,U(:,i));
%         B = vertcat(B,U);
%         C = vertcat(C,C0);
%     end
%     A0 = A0(2:end);
%     A = repmat(A0,1,r);
%     N(:,:,k) = -(A.*B)'*C;
% end
