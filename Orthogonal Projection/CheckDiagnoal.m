function f=CheckDiagnoal(U)
%CHECHDIAGONAL returns 1 and 0 to show if the input U is a diagonal matrix
n=size(U,2);
for i=1:n
    for j=1:n
        if i~=j && norm(U(i,j))>10^(-15)
            f=0;
            return
        end
    end
end
f=1;