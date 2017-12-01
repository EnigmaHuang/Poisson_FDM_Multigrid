function A = Poisson3D_7pt_GenMat(p)
% Form the coefficient matrix for 3D Poisson equation with 7-point difference
% Each edge has 2^p-1 inner grid points and 2 boundary points
	n  = 2^p - 1;      % Number of inner grid points on an edge
	N  = n * n * n;    % Number of inner grid points on (0,1)*(0,1)*(0,1)
	n2 = n * n;
	
	A = sparse(N, N);
	for kz = 1 : n
		for jy = 1 : n
			for ix = 1 : n
				k = (kz - 1) * n2 + (jy - 1) * n + ix;
				A(k, k) = 6;
				if (ix < n) A(k, k + 1)  = -1; end
				if (ix > 1) A(k, k - 1)  = -1; end
				if (jy < n) A(k, k + n)  = -1; end
				if (jy > 1) A(k, k - n)  = -1; end
				if (kz < n) A(k, k + n2) = -1; end
				if (kz > 1) A(k, k - n2) = -1; end
			end
		end
	end
end