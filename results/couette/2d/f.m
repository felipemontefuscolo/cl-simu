function ff = f(A)
n = size(A,1);
m = size(A,2);

ff = zeros(n-1,m);

for	i=1:n-1
	for j=1:m
		ff(i,j) = log2(A(i+1,j)/A(i,j));
	end
end
