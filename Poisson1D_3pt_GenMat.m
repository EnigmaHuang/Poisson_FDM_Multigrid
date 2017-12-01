function A = Laplace1D_GenMat(n)
	A = sparse(n, n);
	
	A(1, 1) =  2.0;
	A(1, 2) = -1.0;
	for i = 2 : n - 1
		A(i, i - 1) = -1.0;
		A(i, i)     =  2.0;
		A(i, i + 1) = -1.0;
	end
	A(n, n)     =  2.0;
	A(n, n - 1) = -1.0;
end