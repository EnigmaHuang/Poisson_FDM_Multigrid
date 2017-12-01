function A = Poisson2D_5pt_GenMat(p)
% Form the coefficient matrix for 2D Poisson equation with 5-point difference
% Each edge has 2^p-1 inner grid points and 2 boundary points
	n = 2^p - 1;      % Number of inner grid points on an edge
	N = n * n;        % Number of inner grid points on (0,1)*(0,1)
	
	A = sparse(N, N); 
	for k = 1 : N
		ix = mod((k - 1), n)    + 1;
		jy = floor((k - 1) / n) + 1;
		
		A(k, k) = 4;
		if (ix < n) A(k, k + 1) = -1; end
		if (ix > 1) A(k, k - 1) = -1; end
		if (jy < n) A(k, k + n) = -1; end
		if (jy > 1) A(k, k - n) = -1; end
	end
end