function [REG, f] = bfPDCCHinterleaver(N_REG_CORESET, L, R, nshift, varargin)
% п.7.3.2.2 TS 38.211
% функция bfPDCCHinterleaver возвращает индексы f для CCE 
% и индексы REG(resource-element groups) 
%(CCE состоят из 6
% REG(resource-element groups) каждый из которых состоит из 12 ресурс 
% элементов (ресурс элемент - одна ячейка на сетке ресурсов)
% N_REG_CORESET количество ресурс блоков для текущего CORESET
% копирование CCE в REG может выполняться как перемешиванием, так и без
% interleaver = 1 активирует функцию перемешивания, перемешивание осуществляется
% кусками REG bundle L - урвень агрегации, параметр верхнего уровня, 
% задает длину REG bundle это число определяет сколько подряд идущих REG
% будет использоваться для копирвания CCE 
% R размер интерливера R = {0,3,6}, задается парамером верхнего уровня
% nshift = N_ID_cell при конфигурации CORESET для PBCH или SIB1
% nshift?{0,1,...,274} в ином случае

% REG представляет собой двумерный массив размерности [6, N_REG_CORESET/6]
% и выдает индексы REG для каждого элемента CCE, всего N_REG_CORESET/6 CСE 
REG = zeros(6, N_REG_CORESET/6);

narginchk(4,5);
if nargin==4
   interleaver = 0;
else
   interleaver = varargin{1};
end

if interleaver 
  % копирования осуществляется с замешиванием  
  C = N_REG_CORESET/(L*R);
  f = zeros(C*R, 1);
  for r = 0 : R - 1
      for c = 0 : C - 1
         j = c*R + r + 1;
         f(j,1) = mod((r*C + c + nshift),(N_REG_CORESET/L));
      end
  end
  
  for j = 0 : N_REG_CORESET/6 - 1
    for k = 0 : 6/L - 1
        for l = 0 : L - 1
            REG(k*L + l + 1, j + 1) = f(6*j/L + k + 1)*6 + l;
        end
    end
  end
else
  % копирования осуществляется без замешивания    
  f = zeros(N_REG_CORESET/6, 1);
  for j = 0 : N_REG_CORESET/6 - 1
     f(j + 1, 1) = j;
     for l = 0 : 5
         REG(l + 1, j + 1) = j * 6 + l;
     end
  end 
end

end