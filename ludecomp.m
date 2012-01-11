% Decomposicao LU sem pivoteamento
% Autor: Afonso Paiva
% Data: 07/02/2010
% Verao 2010 -- ICMC-USP
% Input: Matriz quadrada A(nxn).
% Output: Matrizes triangulares inferior L e superior U de A = LU.

function [L,U]=ludecomp(A)
n=size(A,1);
L=eye(n); U=zeros(n);

for k=1:n
    for j=k:n
        U(k,j)= A(k,j)-L(k,1:k-1)*U(1:k-1,j);
    end
    for i=k+1:n
        L(i,k)=(A(i,k)-L(i,1:k-1)*U(1:k-1,k))/U(k,k);
    end
end

