function out = bfGetKron(N)
n = log2(N);
if n ~= round(n)
    error('N should be a power of 2');
end
  G2 = [1 0; 1 1];
  G_N = 1;
  for i = 1 : n
    G_N = kroneker(G_N,G2);
  end 
  out = G_N;
end

function out = kroneker(A,B)
% ссылка на страницу https://wiki2.org/ru/%D0%9F%D1%80%D0%BE%D0%B8%D0%B7%D0%B2%D0%B5%D0%B4%D0%B5%D0%BD%D0%B8%D0%B5_%D0%9A%D1%80%D0%BE%D0%BD%D0%B5%D0%BA%D0%B5%D1%80%D0%B0
% Если A - матрица размера  m?n, BB — матрица размера p?q, тогда 
% произведение Кронекера есть блочная матрица размера mp?nq

  [m,n] = size(A);
  [p,q] = size(B);
  res = zeros(m*p,n*q);
  
  for mm = 0 : m - 1
    for nn = 0 : n - 1
      for pp = 0 : p - 1
        for qq = 0 : q - 1
            res(mm*p + pp + 1, (nn*q) + qq + 1) = A(mm + 1,nn + 1) * B(pp + 1, qq + 1) ;
        end 
      end
    end
  end
  
  out = res;
end