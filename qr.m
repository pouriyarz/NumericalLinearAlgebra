function [Q, R] = qr(A)
  [m, n] = size(A);

  if (m > n)
    R = zeros(m, n);
    Q = eye(m, m);
    z = zeros(m);

    % Main Program for QR Algorithm, m > n
    for j = 1:n
        % Reflection of column vector of A
        y = A(j:end,j);
        sgn = sign(A(j,j));
        if sgn == 0
          sgn = 1;
        end
        w = y + (sgn * norm(y) * eye(size(y,1),1));
        v = w / norm(w);
        Dummy = 2 * (v * v');
        z(j:end,j:end) = Dummy;
        % Generating Househoulder Matrix
        H = eye(m) - z;
        % Calculating new matrix A using H*A
        A = H * A;
        % Calculating orthogonal matrix, Q using Q=H1*H2*....*Hn
        Q = Q * H;
        z = zeros(m);
    end
  
  else
    % Main Program for QR Algorithm, m < n 
    R = zeros(m, n);
    Q = eye(m, m);
    z = zeros(m);
    for j = 1:m-1
      % Reflection of column vector of A
      y = A(j:end,j);
      sgn = sign(A(j,j));
      if sgn == 0
        sgn = 1;
      end
      w = y +sgn*norm(y)*eye(size(y,1),1);
      v = w/norm(w);
      Dummy = 2*(v*v');
      z(j:end,j:end) = Dummy;
      % Generating Househoulder Matrix
      H = eye(m) - z;
      % Calculating new matrix A using H*A
      A = H*A;
      % Calculating orthogonal matrix, Q using Q=H1*H2*....*Hn
      Q = Q * H;
      z = zeros(m);
    end
  end
  % Forming the R matrix, R = A
  for i = 1:m
    for j = i:n
      R(i,j)=A(i,j);
    end
  end
  % QR decomposition of A
  fprintf('\n  QR DECOMPOSITION OF A :\n');
end
