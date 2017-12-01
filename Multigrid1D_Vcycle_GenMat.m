function [A_list, R_list, max_level] = Multigrid1D_Vcycle_GenMat(A, direct_n)
	n = size(A, 1);
	
	A_list = {};
	R_list = {};
	
	level = 1;
	A_list(level) = {A};
	
	while (n > direct_n)
		coarse_n = floor((n - 1) / 2);
		
		% Construct full-weighted restriction operator 
		R = sparse(coarse_n, n);
		for i = 1 : coarse_n
			col = 2 * i - 1;
			R(i, col)     = 0.25;
			R(i, col + 1) = 0.5;
			R(i, col + 2) = 0.25;
		end
		
		R_list(level) = {R};
		
		P = 2 * R';
		A = R * A * P;
		
		level = level + 1;
		A_list(level) = {A};
		n = coarse_n;
	end
	
	max_level = level;
end