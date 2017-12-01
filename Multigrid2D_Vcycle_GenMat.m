function [A_list, R_list, max_level] = Multigrid2D_Vcycle_GenMat(A, direct_N)
	N = size(A, 1);
	n = floor(sqrt(N));
	coarse_n = floor((n - 1) / 2);
	coarse_N = coarse_n * coarse_n;
	
	A_list = {};
	R_list = {};
	
	level = 1;
	A_list(level) = {A};
	
	while (N > direct_N)
		coarse_n = floor((n - 1) / 2);
		coarse_N = coarse_n * coarse_n;
		
		% Construct restriction operators
		R = sparse(coarse_N, N);
		k = 0;
		for jy = 2 : 2 : n
			for ix = 2 : 2 : n
				k = k + 1;
				fine_grid_k = (jy - 1) * n + ix;
				
				R(k, fine_grid_k - n - 1) = 0.0625;
				R(k, fine_grid_k - n   )  = 0.125;
				R(k, fine_grid_k - n + 1) = 0.0625;
				
				R(k, fine_grid_k - 1)     = 0.125;
				R(k, fine_grid_k    )     = 0.25;
				R(k, fine_grid_k + 1)     = 0.125;
				
				R(k, fine_grid_k + n - 1) = 0.0625;
				R(k, fine_grid_k + n   )  = 0.125;
				R(k, fine_grid_k + n + 1) = 0.0625;
			end
		end
		
		R_list(level) = {R};
		
		P = 4 * R';
		A = R * A * P;
		
		level = level + 1;
		A_list(level) = {A};
		n = coarse_n;
		N = coarse_N;
	end
	max_level = level;
end