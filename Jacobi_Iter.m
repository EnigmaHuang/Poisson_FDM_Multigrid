function [x, iter_cnt] = Jacobi_Iter(A, b, iter_eps, max_iter, x0)
	if size(A, 1) ~= size(A, 2)
        x = -1;
        iter_cnt = 0;
        return;
    end
	
	n = size(A, 1);
	iter_cnt = 0;
	
	if (nargin < 3)	iter_eps = 1e-10; end
	if (nargin < 4)	max_iter = 1000;  end
    if (nargin < 5) 
		x = zeros(n, 1);
	else
		x = x0;
	end
	
	L = -tril(A, -1);
	U = -triu(A, 1);
	D = diag(diag(A));
    B = inv(D) * (L + U);
	g = inv(D) * b;
	
	err = 1;
	while ((err >= iter_eps) && (iter_cnt < max_iter))
		x0 = x;
		x = B * x0 + g;
		iter_cnt = iter_cnt + 1;
		err = max(abs(x0 - x));
	end
end