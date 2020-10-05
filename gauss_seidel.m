function [x] = gauss_seidel(A, b, x, k)
  [m,n] = size(A);
  if m ~= n
    disp('Matris Morabai Nist');
    return;
  end
  x = x';
  new_x = x;
  D = zeros(n);
  D = diag(diag(A));
  L = tril(A) - D;
  U = triu(A) - D;
  Bgs = -inv(D + L) * U;
  Cgs = inv(D + L) * b';
  for i=1:k
    new_x = (Bgs * x) + Cgs;
    x = new_x;
  end
end
