function [A_list, R_list, max_level] = Multigrid3D_Vcycle_GenMat(A, direct_N)
	N = size(A, 1);
	n = round(N^(1 / 3));
	
	A_list = {};
	R_list = {};
	
	level = 1;
	A_list(level) = {A};
	
	while (N > direct_N)
		coarse_n = floor((n - 1) / 2);
		coarse_N = coarse_n * coarse_n * coarse_n;
		
		% Construct restriction operators
		R = sparse(coarse_N, N);
		k = 0;
		n2 = n * n;
		for kz = 2 : 2 : n
			for jy = 2 : 2 : n
				for ix = 2 : 2 : n
					k = k + 1;
					fine_grid_k = (kz - 1) * n2 + (jy - 1) * n + ix;
					
					% Plane z = kz
					R(k, fine_grid_k - n - 1) = 1 / 32;
					R(k, fine_grid_k - n   )  = 1 / 16;
					R(k, fine_grid_k - n + 1) = 1 / 32;
					
					R(k, fine_grid_k - 1)     = 1 / 16;
					R(k, fine_grid_k    )     = 1 / 8;
					R(k, fine_grid_k + 1)     = 1 / 16;
					
					R(k, fine_grid_k + n - 1) = 1 / 32;
					R(k, fine_grid_k + n   )  = 1 / 16;
					R(k, fine_grid_k + n + 1) = 1 / 32;
					
					% Plane z = kz - 1
					R(k, fine_grid_k - n2 - n - 1) = 1 / 64;
					R(k, fine_grid_k - n2 - n   )  = 1 / 32;
					R(k, fine_grid_k - n2 - n + 1) = 1 / 64;
					
					R(k, fine_grid_k - n2 - 1)     = 1 / 32;
					R(k, fine_grid_k - n2    )     = 1 / 16;
					R(k, fine_grid_k - n2 + 1)     = 1 / 32;
					
					R(k, fine_grid_k - n2 + n - 1) = 1 / 64;
					R(k, fine_grid_k - n2 + n   )  = 1 / 32;
					R(k, fine_grid_k - n2 + n + 1) = 1 / 64;
					
					% Plane z = kz + 1
					R(k, fine_grid_k + n2 - n - 1) = 1 / 64;
					R(k, fine_grid_k + n2 - n   )  = 1 / 32;
					R(k, fine_grid_k + n2 - n + 1) = 1 / 64;
					
					R(k, fine_grid_k + n2 - 1)     = 1 / 32;
					R(k, fine_grid_k + n2    )     = 1 / 16;
					R(k, fine_grid_k + n2 + 1)     = 1 / 32;
					
					R(k, fine_grid_k + n2 + n - 1) = 1 / 64;
					R(k, fine_grid_k + n2 + n   )  = 1 / 32;
					R(k, fine_grid_k + n2 + n + 1) = 1 / 64;
				end
			end
		end
		
		R_list(level) = {R};
		
		P = 8 * R';
		A = R * A * P;
		
		level = level + 1;
		A_list(level) = {A};
		n = coarse_n;
		N = coarse_N;
	end
	max_level = level;
end