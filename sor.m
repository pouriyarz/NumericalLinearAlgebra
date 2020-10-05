function [x] = sor(A, b, x, k)
  [m, n] = size(A);
  if m ~= n
    disp('Matris Morabai Nist');
    return;
  end
  x = x';
  b = b';
  new_x = x;
  D = zeros(n);
  D = diag(diag(A));
  L = tril(A) - D;
  U = triu(A) - D;
  
  Bj = inv(D) * (L + U);
  Rj = max(abs(eig(Bj)));
  Rgs = Rj^2;
  if(Rgs > 1)
    disp('ERROR: Hamgera Nist!');
    return;
  end
  w = 2 / (1 + sqrt(1 - Rgs));
  Bs=(inv(D + (w * L))) * (((1 - w) * D) -(w * U));
  Cs=(inv(D + (w * L))) * w * b;
  for i=1:k
    new_x = (Bs * x) + Cs;
    x = new_x ;
  end
end
