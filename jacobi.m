function [x] = jacobi(A, b, x, k)
  [m,n] = size(A);
  if m ~= n
    disp('Matris Morabai Nist');
    return;
  end
  x = x';     # x ra bayad taranahade konim
  new_x = x;  # x jadid baraye estefade dar algorithm
  D = zeros(n);
  D = diag(diag(A));
  L = tril(A) - D;
  U = triu(A) - D;
  Bj = -inv(D) * (L + U);
  Cj = inv(D) * b';
  for i=1:k
    new_x = (Bj * x) + Cj;
    x = new_x ;
  end
end
