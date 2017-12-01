function test_Poisson2D_5pt(p)
	n = 2^p - 1;      % Number of inner grid points on an edge
	N = n * n;        % Number of inner grid points on (0,1)*(0,1)
	rng(n);
	fprintf('Using %d * %d square initial grid\n', n, n);
	A = Poisson2D_5pt_GenMat(p);
	b = rand(N, 1) - 0.5;
	[x, vc_cnt] = Multigrid_Solver(A, b, 2);
end