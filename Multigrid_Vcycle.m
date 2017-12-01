function x = Multigrid_Vcycle(level, A_list, R_list, b, x0, direct_n, PR_coef, smoother, pre_steps, pos_steps)
% Multigrid V-clcye
% level     : The current level (initial is 1)
% A_list    : The array of coefficient matrices on each level
% R_list    : The array of restriction operators on each level
% b         : The right hand side
% x0        : Initial guess
% direct_n  : Threshold for directly solving A_list(level) * x = b
% PR_coef   : The coefficient constant between restriction and prolongation
% smoother  : Function handle for a iterative method as a smoother
% pre_steps : Number of iterations in the pre-smoothing
% pos_steps : Number of iterations in the post-smoothing
	if (nargin < 8)  smoother  = @GS_Iter; end
	if (nargin < 9)  pre_steps = 1;        end
	if (nargin < 10) pos_steps = 1;        end
	
	% Load coefficient matrix
	A = cell2mat(A_list(level));
	
	% If the problem is small enough, solve it directly
	n = size(b, 1);
	if (n <= direct_n)
		x = A \ b;
		return;
	end

	% Pre-smoothing
	x = smoother(A, b, 1e-14, pre_steps, x0);
	
	% Load restriction operator and construct interpolation operator
	R = cell2mat(R_list(level));
	P = R' * PR_coef;
	coarse_n = size(R, 1);
	
	% Compute residual and transfer to coarse grid
	r   = b - A * x;
	r_H = R * r;
	
	% Solve coarse grid problem recursively
	x0  = zeros(coarse_n, 1);
	e_H = Multigrid_Vcycle(level + 1, A_list, R_list, r_H, x0, direct_n, PR_coef, smoother, pre_steps, pos_steps);
	
	% Transfer error to fine grid and correct
	x = x + P * e_H;
	
	% Post-smoothing
	x = smoother(A, b, 1e-14, pos_steps, x);
end