function [x, vcycle_cnt, res_norm] = Multigrid_Solver(A, b, dim, smoother, pre_steps, pos_steps, rn_tol)
% Multigrid solver for A * x = b on n-dimension cude grid
% A         : The inital coefficient matrix
% b         : The right hand side
% dim       : The dimension of the grid (1, 2 or 3)
% smoother  : Function handle for a iterative method as a smoother
% pre_steps : Number of iterations in the pre-smoothing
% pos_steps : Number of iterations in the post-smoothing
% rn_tol    : The tolerance of the relative residual norm
	if (nargin < 4) smoother  = @GS_Iter; end
	if (nargin < 5) pre_steps = 1;	      end
	if (nargin < 6) pos_steps = 1;	      end
	if (nargin < 7) rn_tol    = 1e-10;    end
	
	n  = size(A, 1);
	x  = zeros(n, 1);
	rn = norm(b);
	vcycle_cnt  = 0;
	res_norm(1) = rn;
	rn_stop     = rn * rn_tol;
	
	if (dim == 1) 
		PR_coef  = 2;
		direct_n = 16; 
	end
	if (dim == 2) 
		PR_coef  = 4;
		direct_n = 7 * 7; 
	end
	if (dim == 3) 
		PR_coef  = 8;
		direct_n = 3 * 3 * 3; 
	end
	
	% Generate coefficient matrices and restriction operators of each level at once
	tic;
	if (dim == 1)
		[A_list, R_list, max_level] = Multigrid1D_Vcycle_GenMat(A, direct_n);
	end
	if (dim == 2)
		[A_list, R_list, max_level] = Multigrid2D_Vcycle_GenMat(A, direct_n);
	end
	if (dim == 3)
		[A_list, R_list, max_level] = Multigrid3D_Vcycle_GenMat(A, direct_n);
	end
	gm_t = toc;
	% Repeat V-cycle until converge
	tic;
	while (rn > rn_stop)
		x  = Multigrid_Vcycle(1, A_list, R_list, b, x, direct_n, PR_coef, smoother, pre_steps, pos_steps);
		r  = b - A * x;
		rn = norm(r, 2);
		vcycle_cnt = vcycle_cnt + 1;
		res_norm(vcycle_cnt + 1) = rn;
	end
	vcyc_t = toc;
	
	fprintf('Matrices generating wall-time = %f (s)\n', gm_t);
	fprintf('V-cycle solver      wall-time = %f (s)\n', vcyc_t);
	fprintf('Performed V-cycles = %d\n', vcycle_cnt);
	fprintf('||b - A * x||_2    = %e\n', rn);
end