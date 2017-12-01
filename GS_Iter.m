function [x, iter_cnt] = GS_Iter(A, b, iter_eps, max_iter, x0)
	if size(A, 1) ~= size(A, 2)
		x = -1;
		iter_cnt = 0;
		return;
	end
	
	n = size(A, 1);
	iter_cnt = 0;
	err = 1;
	
	if (nargin < 3)	iter_eps = 1e-10; end
	if (nargin < 4)	max_iter = 1000;  end
	if (nargin < 5) 
		x = zeros(n, 1) + iter_eps;
	else
		x = x0;
	end

	while ((err >= iter_eps) && (iter_cnt < max_iter))
		x0 = x;
		
		for i = 1 : n
			A_ii = A(i, i);
			x(i) = b(i) - A(i, :) * x + A_ii * x(i);
			x(i) = x(i) / A_ii;
		end
		
		iter_cnt = iter_cnt + 1;
		err = max(abs(x0 - x));
	end
end