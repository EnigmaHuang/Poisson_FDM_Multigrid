function test_Poisson3D(p)
	n = 2^p - 1;      % Number of inner grid points on an edge
	N = n * n * n;    % Number of inner grid points on (0,1)*(0,1)*(0,1)
	rng(n);
	fprintf('Using %d * %d * %d initial cube\n', n, n, n);
	A = Poisson3D_7pt_GenMat(p);
	b = rand(N, 1) - 0.5;
	[x, vc_cnt] = Multigrid_Solver(A, b, 3);
end