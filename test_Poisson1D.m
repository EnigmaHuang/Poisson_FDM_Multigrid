function test_Poisson1D(p)
	n = 2^p - 1;
	rng(n);
	fprintf('Using %d initial grid points\n', n);
	A = Poisson1D_3pt_GenMat(n);
	b = rand(n, 1) - 0.5;
	[x, vc_cnt] = Multigrid_Solver(A, b, 1);
end